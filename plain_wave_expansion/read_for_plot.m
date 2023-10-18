a=load('wk--121NiAl.mat',"-ascii")
figure;
for i=1:121
    c=a(:,i);
    j=1:60;
    plot(j,c,'r.');
    hold on;
end
xlabel('wave vector')
ylabel('normalized frequency')
axis([1,60,0,1]);
