function L=compute_length(NN,G)

L = sqrt((NN(1,G(2,:))-NN(1,G(1,:))).^2+...
         (NN(2,G(2,:))-NN(2,G(1,:))).^2+...
         (NN(3,G(2,:))-NN(3,G(1,:))).^2);

end