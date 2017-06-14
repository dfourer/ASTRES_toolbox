function theta2=convert_angle(theta1)

theta2=theta1;
ind1=find(theta1>0);
theta2(ind1)=pi/2-theta1(ind1);

ind2=find(theta1<0);
theta2(ind2)=-pi/2-theta1(ind2);