function snr_value = compute_snr(original, recovered)
    noise = original - recovered;
    snr_value = 10 * log10(sum(original(:).^2) / sum(noise(:).^2));
end