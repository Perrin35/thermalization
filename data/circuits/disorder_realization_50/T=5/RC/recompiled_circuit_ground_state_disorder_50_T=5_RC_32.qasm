OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.1640373) q[0];
sx q[0];
rz(-1.9174175) q[0];
sx q[0];
rz(1.8608215) q[0];
rz(2.3692853) q[1];
sx q[1];
rz(-2.5025855) q[1];
sx q[1];
rz(-0.26689902) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.229202) q[0];
sx q[0];
rz(-1.664219) q[0];
sx q[0];
rz(-1.6673198) q[0];
x q[1];
rz(2.2592531) q[2];
sx q[2];
rz(-1.7710476) q[2];
sx q[2];
rz(0.74742521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3223826) q[1];
sx q[1];
rz(-1.1448082) q[1];
sx q[1];
rz(-2.3374578) q[1];
rz(1.4679174) q[3];
sx q[3];
rz(-1.130338) q[3];
sx q[3];
rz(-0.48635024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.018517) q[2];
sx q[2];
rz(-0.96394959) q[2];
sx q[2];
rz(-0.71259552) q[2];
rz(2.9739042) q[3];
sx q[3];
rz(-2.2919787) q[3];
sx q[3];
rz(2.6970421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0406822) q[0];
sx q[0];
rz(-1.7829144) q[0];
sx q[0];
rz(1.6810625) q[0];
rz(-1.4617317) q[1];
sx q[1];
rz(-2.0748383) q[1];
sx q[1];
rz(2.0251958) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76976639) q[0];
sx q[0];
rz(-0.049749181) q[0];
sx q[0];
rz(-1.0509963) q[0];
rz(-pi) q[1];
rz(-0.15123151) q[2];
sx q[2];
rz(-2.1215212) q[2];
sx q[2];
rz(-3.0557003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9560686) q[1];
sx q[1];
rz(-1.9466234) q[1];
sx q[1];
rz(1.2898499) q[1];
rz(0.41023572) q[3];
sx q[3];
rz(-1.6950775) q[3];
sx q[3];
rz(-1.6491485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1560893) q[2];
sx q[2];
rz(-2.6458793) q[2];
sx q[2];
rz(-0.25225857) q[2];
rz(-1.0559399) q[3];
sx q[3];
rz(-0.85774937) q[3];
sx q[3];
rz(2.1411538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58241874) q[0];
sx q[0];
rz(-2.2370339) q[0];
sx q[0];
rz(-2.3140123) q[0];
rz(-2.6170392) q[1];
sx q[1];
rz(-1.8962212) q[1];
sx q[1];
rz(0.96447271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37999718) q[0];
sx q[0];
rz(-0.58172136) q[0];
sx q[0];
rz(-1.6935478) q[0];
x q[1];
rz(3.0542637) q[2];
sx q[2];
rz(-1.9448408) q[2];
sx q[2];
rz(2.5279074) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71389329) q[1];
sx q[1];
rz(-1.5212219) q[1];
sx q[1];
rz(2.5109641) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6286704) q[3];
sx q[3];
rz(-2.0116968) q[3];
sx q[3];
rz(2.1564275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7437637) q[2];
sx q[2];
rz(-0.22481329) q[2];
sx q[2];
rz(-1.1986097) q[2];
rz(-1.648929) q[3];
sx q[3];
rz(-2.1677833) q[3];
sx q[3];
rz(0.57884136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32597932) q[0];
sx q[0];
rz(-0.15376832) q[0];
sx q[0];
rz(1.9258668) q[0];
rz(0.99023306) q[1];
sx q[1];
rz(-2.4001887) q[1];
sx q[1];
rz(-2.5726817) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7051681) q[0];
sx q[0];
rz(-1.8537844) q[0];
sx q[0];
rz(-2.1678336) q[0];
x q[1];
rz(0.26428916) q[2];
sx q[2];
rz(-1.2131872) q[2];
sx q[2];
rz(2.8823095) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.68142707) q[1];
sx q[1];
rz(-2.6003464) q[1];
sx q[1];
rz(0.98249225) q[1];
x q[2];
rz(-0.91671555) q[3];
sx q[3];
rz(-0.93610969) q[3];
sx q[3];
rz(-0.095528729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27900055) q[2];
sx q[2];
rz(-2.180884) q[2];
sx q[2];
rz(-0.29435364) q[2];
rz(-0.66458464) q[3];
sx q[3];
rz(-2.3817101) q[3];
sx q[3];
rz(-1.6871281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6890474) q[0];
sx q[0];
rz(-1.4280467) q[0];
sx q[0];
rz(0.31387615) q[0];
rz(0.9017871) q[1];
sx q[1];
rz(-1.7867463) q[1];
sx q[1];
rz(-1.6759759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2814502) q[0];
sx q[0];
rz(-1.6897795) q[0];
sx q[0];
rz(-0.02597474) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1526665) q[2];
sx q[2];
rz(-2.8684432) q[2];
sx q[2];
rz(1.7371617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5295388) q[1];
sx q[1];
rz(-1.1837237) q[1];
sx q[1];
rz(-1.6516634) q[1];
rz(0.05081049) q[3];
sx q[3];
rz(-0.63119315) q[3];
sx q[3];
rz(2.7630591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9471211) q[2];
sx q[2];
rz(-0.7496382) q[2];
sx q[2];
rz(2.0571902) q[2];
rz(0.32367745) q[3];
sx q[3];
rz(-2.4249228) q[3];
sx q[3];
rz(2.9173541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502458) q[0];
sx q[0];
rz(-2.8960189) q[0];
sx q[0];
rz(-2.7299951) q[0];
rz(-0.72548524) q[1];
sx q[1];
rz(-1.8975703) q[1];
sx q[1];
rz(1.4010319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4759504) q[0];
sx q[0];
rz(-0.67806292) q[0];
sx q[0];
rz(-0.31291385) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13757041) q[2];
sx q[2];
rz(-2.3964632) q[2];
sx q[2];
rz(2.957475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7606719) q[1];
sx q[1];
rz(-0.47539818) q[1];
sx q[1];
rz(0.74785079) q[1];
x q[2];
rz(-1.6468923) q[3];
sx q[3];
rz(-2.020936) q[3];
sx q[3];
rz(-1.744818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1650042) q[2];
sx q[2];
rz(-0.62825957) q[2];
sx q[2];
rz(1.4264301) q[2];
rz(-0.12380883) q[3];
sx q[3];
rz(-2.0721469) q[3];
sx q[3];
rz(2.6482705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2774778) q[0];
sx q[0];
rz(-1.1461698) q[0];
sx q[0];
rz(1.3364828) q[0];
rz(0.9388963) q[1];
sx q[1];
rz(-2.0195596) q[1];
sx q[1];
rz(0.28405651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1079252) q[0];
sx q[0];
rz(-1.9104275) q[0];
sx q[0];
rz(0.64396283) q[0];
rz(-pi) q[1];
rz(-0.16347537) q[2];
sx q[2];
rz(-1.7942085) q[2];
sx q[2];
rz(-2.7976409) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3160243) q[1];
sx q[1];
rz(-3.1085759) q[1];
sx q[1];
rz(1.2570639) q[1];
rz(1.7098996) q[3];
sx q[3];
rz(-1.3918607) q[3];
sx q[3];
rz(-0.034497189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9629918) q[2];
sx q[2];
rz(-2.4602349) q[2];
sx q[2];
rz(-0.94092384) q[2];
rz(-3.1407147) q[3];
sx q[3];
rz(-1.9160756) q[3];
sx q[3];
rz(0.0860478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0451999) q[0];
sx q[0];
rz(-2.2470076) q[0];
sx q[0];
rz(2.9659502) q[0];
rz(1.9215709) q[1];
sx q[1];
rz(-0.86985391) q[1];
sx q[1];
rz(-0.87179914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891465) q[0];
sx q[0];
rz(-1.8348357) q[0];
sx q[0];
rz(-0.15588481) q[0];
rz(-1.9398676) q[2];
sx q[2];
rz(-0.82316527) q[2];
sx q[2];
rz(1.6047503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2928472) q[1];
sx q[1];
rz(-2.1715144) q[1];
sx q[1];
rz(1.377618) q[1];
x q[2];
rz(-0.11856793) q[3];
sx q[3];
rz(-0.55988479) q[3];
sx q[3];
rz(2.9836451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.72622) q[2];
sx q[2];
rz(-1.8612334) q[2];
sx q[2];
rz(-1.0225164) q[2];
rz(2.7546049) q[3];
sx q[3];
rz(-1.756668) q[3];
sx q[3];
rz(-2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33685327) q[0];
sx q[0];
rz(-1.1664718) q[0];
sx q[0];
rz(0.26300305) q[0];
rz(2.8291342) q[1];
sx q[1];
rz(-1.0954906) q[1];
sx q[1];
rz(-1.9591263) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1240886) q[0];
sx q[0];
rz(-0.54709439) q[0];
sx q[0];
rz(1.1897683) q[0];
rz(-0.55652852) q[2];
sx q[2];
rz(-0.91224837) q[2];
sx q[2];
rz(-1.5040894) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14543372) q[1];
sx q[1];
rz(-2.7163032) q[1];
sx q[1];
rz(0.49844679) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1402604) q[3];
sx q[3];
rz(-2.0234442) q[3];
sx q[3];
rz(-2.4695726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.075228127) q[2];
sx q[2];
rz(-1.790739) q[2];
sx q[2];
rz(2.8933375) q[2];
rz(-2.9634641) q[3];
sx q[3];
rz(-2.0592561) q[3];
sx q[3];
rz(2.3586912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7149776) q[0];
sx q[0];
rz(-2.7474032) q[0];
sx q[0];
rz(1.0618427) q[0];
rz(-1.8408403) q[1];
sx q[1];
rz(-1.6164833) q[1];
sx q[1];
rz(-1.1311857) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6908195) q[0];
sx q[0];
rz(-1.1191219) q[0];
sx q[0];
rz(-1.0502688) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50297351) q[2];
sx q[2];
rz(-1.9679234) q[2];
sx q[2];
rz(3.0850038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0319236) q[1];
sx q[1];
rz(-1.1732035) q[1];
sx q[1];
rz(-1.470674) q[1];
rz(1.3374109) q[3];
sx q[3];
rz(-2.1595229) q[3];
sx q[3];
rz(-1.0860541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0807557) q[2];
sx q[2];
rz(-2.8751825) q[2];
sx q[2];
rz(-2.5401435) q[2];
rz(-2.1001749) q[3];
sx q[3];
rz(-1.6626549) q[3];
sx q[3];
rz(-1.3781594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9926485) q[0];
sx q[0];
rz(-0.61275488) q[0];
sx q[0];
rz(-2.3339094) q[0];
rz(-3.0638937) q[1];
sx q[1];
rz(-0.54010375) q[1];
sx q[1];
rz(1.7355951) q[1];
rz(-0.97066047) q[2];
sx q[2];
rz(-1.9117336) q[2];
sx q[2];
rz(-2.5993549) q[2];
rz(0.79176767) q[3];
sx q[3];
rz(-2.0590012) q[3];
sx q[3];
rz(1.7760359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
