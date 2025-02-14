OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5975098) q[0];
sx q[0];
rz(-1.37356) q[0];
sx q[0];
rz(-2.1231667) q[0];
rz(2.6234558) q[1];
sx q[1];
rz(-0.75704804) q[1];
sx q[1];
rz(2.5098324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7512832) q[0];
sx q[0];
rz(-1.3645118) q[0];
sx q[0];
rz(0.12796107) q[0];
x q[1];
rz(-1.6502981) q[2];
sx q[2];
rz(-2.8379746) q[2];
sx q[2];
rz(2.7111766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0644016) q[1];
sx q[1];
rz(-0.6997515) q[1];
sx q[1];
rz(-1.1279068) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8381836) q[3];
sx q[3];
rz(-0.33743706) q[3];
sx q[3];
rz(-2.9149559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6112001) q[2];
sx q[2];
rz(-0.35553122) q[2];
sx q[2];
rz(1.4872888) q[2];
rz(0.90855956) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(2.1366185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0205883) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(2.9793136) q[0];
rz(0.24457112) q[1];
sx q[1];
rz(-1.2030315) q[1];
sx q[1];
rz(2.8499106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26839089) q[0];
sx q[0];
rz(-1.8531728) q[0];
sx q[0];
rz(1.3703521) q[0];
rz(0.9082011) q[2];
sx q[2];
rz(-2.78763) q[2];
sx q[2];
rz(2.5733054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32714601) q[1];
sx q[1];
rz(-1.8015878) q[1];
sx q[1];
rz(2.3620741) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41145153) q[3];
sx q[3];
rz(-1.9570159) q[3];
sx q[3];
rz(0.4377818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4633999) q[2];
sx q[2];
rz(-0.86898154) q[2];
sx q[2];
rz(2.0657067) q[2];
rz(1.894527) q[3];
sx q[3];
rz(-1.5970634) q[3];
sx q[3];
rz(-1.6943078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4397044) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(-0.18371789) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-2.3972062) q[1];
sx q[1];
rz(2.9768129) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69075981) q[0];
sx q[0];
rz(-2.1565821) q[0];
sx q[0];
rz(3.1302439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38305958) q[2];
sx q[2];
rz(-0.7128517) q[2];
sx q[2];
rz(1.0733611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91080571) q[1];
sx q[1];
rz(-1.9525098) q[1];
sx q[1];
rz(-0.64818212) q[1];
rz(-2.1042473) q[3];
sx q[3];
rz(-1.8672393) q[3];
sx q[3];
rz(-0.25755778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.9811337) q[2];
sx q[2];
rz(-1.1064233) q[2];
rz(0.70118457) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(0.4655233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3727386) q[0];
sx q[0];
rz(-2.0059858) q[0];
sx q[0];
rz(2.1970774) q[0];
rz(-1.5185897) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(1.6015582) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2006783) q[0];
sx q[0];
rz(-0.8884065) q[0];
sx q[0];
rz(-0.81253482) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11941274) q[2];
sx q[2];
rz(-1.4541199) q[2];
sx q[2];
rz(1.2774955) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9303927) q[1];
sx q[1];
rz(-0.89387776) q[1];
sx q[1];
rz(2.6571214) q[1];
rz(-2.3911658) q[3];
sx q[3];
rz(-1.6035169) q[3];
sx q[3];
rz(-1.0270384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.010178415) q[2];
sx q[2];
rz(-2.5203036) q[2];
sx q[2];
rz(-0.16769257) q[2];
rz(0.0532648) q[3];
sx q[3];
rz(-1.1054509) q[3];
sx q[3];
rz(-1.7648599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57529706) q[0];
sx q[0];
rz(-0.10558858) q[0];
sx q[0];
rz(-0.29933023) q[0];
rz(0.37295595) q[1];
sx q[1];
rz(-1.2495709) q[1];
sx q[1];
rz(1.4704871) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18749593) q[0];
sx q[0];
rz(-2.4135655) q[0];
sx q[0];
rz(-0.67250979) q[0];
x q[1];
rz(0.46385455) q[2];
sx q[2];
rz(-1.4632483) q[2];
sx q[2];
rz(-0.35075089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5467087) q[1];
sx q[1];
rz(-1.0011893) q[1];
sx q[1];
rz(2.1807266) q[1];
rz(-0.94909747) q[3];
sx q[3];
rz(-1.4273564) q[3];
sx q[3];
rz(-2.6196817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2028929) q[2];
sx q[2];
rz(-0.6627658) q[2];
sx q[2];
rz(1.1104442) q[2];
rz(2.2796196) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(1.2601669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2170169) q[0];
sx q[0];
rz(-2.45455) q[0];
sx q[0];
rz(-2.7365015) q[0];
rz(-2.8054667) q[1];
sx q[1];
rz(-1.7205709) q[1];
sx q[1];
rz(-0.95132557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4753805) q[0];
sx q[0];
rz(-2.6586091) q[0];
sx q[0];
rz(0.57421143) q[0];
rz(-pi) q[1];
rz(-1.3589922) q[2];
sx q[2];
rz(-2.2957605) q[2];
sx q[2];
rz(-2.447809) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7860884) q[1];
sx q[1];
rz(-2.5230683) q[1];
sx q[1];
rz(2.7433646) q[1];
x q[2];
rz(1.7634994) q[3];
sx q[3];
rz(-1.2964222) q[3];
sx q[3];
rz(1.5128795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0508017) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(-2.9108099) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(-0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7886605) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(-2.2094862) q[0];
rz(-2.508029) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(-1.8036141) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0377696) q[0];
sx q[0];
rz(-1.4892206) q[0];
sx q[0];
rz(-1.4647116) q[0];
rz(-1.6858844) q[2];
sx q[2];
rz(-1.8427927) q[2];
sx q[2];
rz(-2.4912262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38398582) q[1];
sx q[1];
rz(-1.1771056) q[1];
sx q[1];
rz(-0.91169731) q[1];
x q[2];
rz(-0.058480992) q[3];
sx q[3];
rz(-0.49904682) q[3];
sx q[3];
rz(-2.8935695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6771217) q[2];
sx q[2];
rz(-1.9147583) q[2];
sx q[2];
rz(1.202549) q[2];
rz(0.36457148) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(-0.83474368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4778336) q[0];
sx q[0];
rz(-2.5680225) q[0];
sx q[0];
rz(-2.9649576) q[0];
rz(2.7015576) q[1];
sx q[1];
rz(-1.9562079) q[1];
sx q[1];
rz(2.1766591) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10213908) q[0];
sx q[0];
rz(-2.9055465) q[0];
sx q[0];
rz(-2.5419767) q[0];
x q[1];
rz(-2.1148483) q[2];
sx q[2];
rz(-2.1026582) q[2];
sx q[2];
rz(0.085467664) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31613708) q[1];
sx q[1];
rz(-1.9333922) q[1];
sx q[1];
rz(0.98939244) q[1];
rz(1.8252402) q[3];
sx q[3];
rz(-1.644205) q[3];
sx q[3];
rz(-2.0919098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9289916) q[2];
sx q[2];
rz(-1.5233728) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(0.11659226) q[3];
sx q[3];
rz(-2.8023585) q[3];
sx q[3];
rz(-0.54178437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5803439) q[0];
sx q[0];
rz(-1.2224226) q[0];
sx q[0];
rz(-2.5196581) q[0];
rz(1.7033345) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(2.4780746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2885482) q[0];
sx q[0];
rz(-2.2877734) q[0];
sx q[0];
rz(-1.8621481) q[0];
rz(-pi) q[1];
rz(-1.8542669) q[2];
sx q[2];
rz(-1.6890235) q[2];
sx q[2];
rz(-2.9663939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92163699) q[1];
sx q[1];
rz(-1.8884648) q[1];
sx q[1];
rz(-0.89584535) q[1];
rz(2.8597707) q[3];
sx q[3];
rz(-1.661206) q[3];
sx q[3];
rz(0.91109959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(-0.36250472) q[2];
rz(0.47719657) q[3];
sx q[3];
rz(-1.0537182) q[3];
sx q[3];
rz(-2.0188735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0839194) q[0];
sx q[0];
rz(-2.4102983) q[0];
sx q[0];
rz(-2.2398563) q[0];
rz(-0.67507356) q[1];
sx q[1];
rz(-0.98552862) q[1];
sx q[1];
rz(2.9973082) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57406232) q[0];
sx q[0];
rz(-1.4882548) q[0];
sx q[0];
rz(1.922419) q[0];
x q[1];
rz(0.21964964) q[2];
sx q[2];
rz(-1.4712508) q[2];
sx q[2];
rz(-1.2964028) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9291097) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(-0.65156619) q[1];
x q[2];
rz(1.3809526) q[3];
sx q[3];
rz(-1.1355054) q[3];
sx q[3];
rz(-1.3224885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(-0.20067659) q[2];
rz(1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(2.5698575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5889482) q[0];
sx q[0];
rz(-0.968796) q[0];
sx q[0];
rz(2.0198685) q[0];
rz(-2.6523392) q[1];
sx q[1];
rz(-1.4514634) q[1];
sx q[1];
rz(-1.0101752) q[1];
rz(-1.1560925) q[2];
sx q[2];
rz(-1.6558052) q[2];
sx q[2];
rz(2.2075352) q[2];
rz(-2.6993467) q[3];
sx q[3];
rz(-1.5447692) q[3];
sx q[3];
rz(-1.0039644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
