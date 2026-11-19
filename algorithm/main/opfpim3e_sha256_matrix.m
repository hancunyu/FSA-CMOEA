function hash = opfpim3e_sha256_matrix(X)
    try
        md=java.security.MessageDigest.getInstance('SHA-256');
        sz=uint64(size(X)); md.update(typecast(sz(:),'int8'));
        md.update(typecast(X(:),'int8'));
        d=typecast(md.digest(),'uint8');
        hash=lower(reshape(dec2hex(d,2).',1,[]));
    catch
        hash=sprintf('fallback_%dx%d_%.17g',size(X,1),size(X,2),sum(X(:)));
    end
end





