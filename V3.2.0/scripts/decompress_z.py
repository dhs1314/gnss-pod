"""Decompress Unix .Z (LZW compress) files using correct bit ordering."""
from pathlib import Path


def decompress_z(data):
    """Decompress Unix compress (.Z) format (LZW, MSB-first bit order)."""
    if len(data) < 4:
        return data
    if data[0] != 0x1f or data[1] != 0x9d:
        return data

    flags = data[2]
    maxbits = flags & 0x1F
    if maxbits < 9 or maxbits > 16:
        print(f'  WARNING: maxbits={maxbits}, trying 16')
        maxbits = 16

    maxmaxcode = 1 << maxbits
    CLEAR = 256
    FIRST = 257

    result = bytearray()
    table = {i: bytes([i]) for i in range(256)}
    n_bits = 9
    next_code = FIRST

    class BitReader:
        def __init__(self, data, start_bit):
            self.data = data
            self.pos = start_bit  # bit position, 0-indexed from start

        def read(self, n):
            """Read n bits MSB-first."""
            val = 0
            for _ in range(n):
                byte_idx = self.pos // 8
                bit_idx = 7 - (self.pos % 8)
                if byte_idx < len(self.data):
                    bit = (self.data[byte_idx] >> bit_idx) & 1
                    val = (val << 1) | bit
                self.pos += 1
            return val

    bits = BitReader(data, 24)  # start after 3-byte header

    # Read first code
    code = bits.read(n_bits)
    if code >= 256:
        return bytes(result)

    result.extend(table[code])
    prev = table[code]

    while bits.pos + n_bits <= len(data) * 8:
        code = bits.read(n_bits)

        if code == CLEAR:
            table = {i: bytes([i]) for i in range(256)}
            n_bits = 9
            next_code = FIRST
            code = bits.read(n_bits)
            if code >= 256:
                break
            result.extend(table[code])
            prev = table[code]
            continue

        if code < 256:
            entry = table[code]
            result.extend(entry)
            table[next_code] = prev + entry[:1]
            next_code += 1
            prev = entry
        elif code < next_code:
            entry = table[code]
            result.extend(entry)
            table[next_code] = prev + entry[:1]
            next_code += 1
            prev = entry
        elif code == next_code:
            entry = prev + prev[:1]
            result.extend(entry)
            table[next_code] = entry
            next_code += 1
            prev = entry
        else:
            break

        if next_code >= (1 << n_bits) and n_bits < maxbits:
            n_bits += 1

    return bytes(result)


def main():
    dest = Path('data/CODE/2024')
    for fname in ['P1C12404_RINEX.DCB.Z', 'P1C12404.DCB.Z', 'P1P22404.DCB.Z']:
        path = dest / fname
        if not path.exists():
            print(f'{fname}: not found')
            continue
        with open(str(path), 'rb') as f:
            compressed = f.read()

        flags = compressed[2]
        print(f'{fname}: {len(compressed)} bytes, maxbits={flags & 0x1F}')

        data = decompress_z(compressed)
        out_path = dest / fname.replace('.Z', '')
        with open(str(out_path), 'wb') as f:
            f.write(data)
        print(f'  -> {len(data)} bytes')

        try:
            text = data.decode('ascii', errors='replace')
            for line in text.split('\n')[:6]:
                if line.strip():
                    print(f'  | {line.rstrip()[:100]}')
        except Exception:
            pass


if __name__ == '__main__':
    main()
