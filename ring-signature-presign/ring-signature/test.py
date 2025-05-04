import time

from ring_signature import rtScheme


def main():
    rt = rtScheme()
    message = b"Hello, Ring Signature!"
    print("Message: 'Hello, Ring Signature!'\n");


    # Setup
    print("Setup params...")
    start_ns = time.perf_counter_ns()
    rt.setup()
    keygen_time_ms = (time.perf_counter_ns() - start_ns) / 1_000_000
    print(f"✓ Setup params in {keygen_time_ms:.2f} ms\n")

    # Key Generation
    print("Generating keypair...")
    start_ns = time.perf_counter_ns()
    pub, priv = rt.rs_keygen()
    keygen_time_ms = (time.perf_counter_ns() - start_ns) / 1_000_000
    print(f"✓ Keys generated in {keygen_time_ms:.2f} ms\n")
    # print(f"Public key (first 5): {pub[1][0].coefficients[:5]}\n")

    # Signing
    print("Generating signature...")
    start_ns = time.perf_counter_ns()
    try:
        z, c, com = rt.rs_sign(message)
        sign_time_ms = (time.perf_counter_ns() - start_ns) / 1_000_000
        print(f"✓ Signature generated in {sign_time_ms:.2f} ms\n")

        # print("\nSignature Details:")
        # print("-----------------")
        # print(f"Challenge (c): {c.coefficients[:5]}")
        # print(f"Response (z): {z[0].coefficients[:5]}")
        # print(f"Witness (w): {com[0].coefficients[:5]}\n")

        # Verification
        print("Verifying signature...")
        start_ns = time.perf_counter_ns()
        valid = rt.rs_verify(message, (z, c, com))
        verify_time_ms = (time.perf_counter_ns() - start_ns) / 1_000_000

        if valid:
            print("✓ Signature verification successful\n")
        else:
            print("✗ Signature verification failed\n")

        # NIZK Setup
        # start_ns = time.perf_counter_ns()
        # rt.NIZK_Setup()
        # NIZKSP_time_ms = (time.perf_counter_ns() - start_ns) / 1_000_000
        # print(f"NIZK Setup Time: {NIZKSP_time_ms:>8.2f} ms\n")

        print("Generating poof...")
        start_ns = time.perf_counter_ns()
        rt.NIZK_Prove(message)
        prove_time_ms = (time.perf_counter_ns() - start_ns) / 1_000_000
        print(f"NIZK Proving Time: {prove_time_ms:>8.2f} ms\n")
    
        print("Performance Summary:")
        print("------------------")
        print(f"Setup params in {keygen_time_ms:.2f} ms")
        print(f"Key Generation: {keygen_time_ms:>8.2f} ms")
        print(f"Signing Time:   {sign_time_ms:>8.2f} ms")
        # print(f"Verify Time:    {verify_time_ms:>8.2f} ms")
        print(f"NIZK Proving Time: {prove_time_ms:>8.2f} ms\n")
        print(
            f"Total Time:     {(sign_time_ms + prove_time_ms):>8.2f} ms"
        )

    except Exception as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main()
