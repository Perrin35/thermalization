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
rz(2.0848701) q[0];
sx q[0];
rz(4.4216006) q[0];
sx q[0];
rz(13.292424) q[0];
rz(0.59745204) q[1];
sx q[1];
rz(-2.6260881) q[1];
sx q[1];
rz(-2.1093624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7643578) q[0];
sx q[0];
rz(-0.66417686) q[0];
sx q[0];
rz(-0.16420495) q[0];
rz(-2.0886253) q[2];
sx q[2];
rz(-1.6047239) q[2];
sx q[2];
rz(-2.9133493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91019708) q[1];
sx q[1];
rz(-0.33581844) q[1];
sx q[1];
rz(0.37792716) q[1];
rz(0.49555669) q[3];
sx q[3];
rz(-2.2751788) q[3];
sx q[3];
rz(0.92951194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44077474) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(3.0536998) q[2];
rz(1.6739738) q[3];
sx q[3];
rz(-1.2940977) q[3];
sx q[3];
rz(0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217594) q[0];
sx q[0];
rz(-1.8353945) q[0];
sx q[0];
rz(1.4922967) q[0];
rz(-0.78798931) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(0.96510395) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3639587) q[0];
sx q[0];
rz(-2.2978362) q[0];
sx q[0];
rz(1.4476089) q[0];
rz(-pi) q[1];
rz(-2.4048058) q[2];
sx q[2];
rz(-1.4649142) q[2];
sx q[2];
rz(0.69241619) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3841159) q[1];
sx q[1];
rz(-2.6429477) q[1];
sx q[1];
rz(-1.9668786) q[1];
rz(-pi) q[2];
rz(-1.4802259) q[3];
sx q[3];
rz(-1.1027059) q[3];
sx q[3];
rz(0.15453574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1284156) q[2];
sx q[2];
rz(-0.33856496) q[2];
sx q[2];
rz(-1.817912) q[2];
rz(-0.88661083) q[3];
sx q[3];
rz(-1.1611791) q[3];
sx q[3];
rz(2.9429341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.897268) q[0];
sx q[0];
rz(-0.97049814) q[0];
sx q[0];
rz(0.6955198) q[0];
rz(0.8356525) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(2.4998891) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7135237) q[0];
sx q[0];
rz(-1.4260573) q[0];
sx q[0];
rz(0.83065303) q[0];
rz(-pi) q[1];
rz(-2.4303248) q[2];
sx q[2];
rz(-1.4431653) q[2];
sx q[2];
rz(0.39943275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1759626) q[1];
sx q[1];
rz(-1.8847191) q[1];
sx q[1];
rz(2.0759275) q[1];
rz(-pi) q[2];
rz(-2.1459177) q[3];
sx q[3];
rz(-2.174587) q[3];
sx q[3];
rz(-0.010494516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.16331638) q[2];
sx q[2];
rz(-2.0127608) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(1.4000019) q[3];
sx q[3];
rz(-1.497523) q[3];
sx q[3];
rz(-1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9982933) q[0];
sx q[0];
rz(-0.87109733) q[0];
sx q[0];
rz(0.77348462) q[0];
rz(3.0553014) q[1];
sx q[1];
rz(-2.6373865) q[1];
sx q[1];
rz(2.6533244) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9241656) q[0];
sx q[0];
rz(-1.1091976) q[0];
sx q[0];
rz(-1.4946926) q[0];
x q[1];
rz(2.4441798) q[2];
sx q[2];
rz(-2.0048755) q[2];
sx q[2];
rz(3.0212439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1459256) q[1];
sx q[1];
rz(-2.2856376) q[1];
sx q[1];
rz(0.30024528) q[1];
x q[2];
rz(-1.9776439) q[3];
sx q[3];
rz(-1.196612) q[3];
sx q[3];
rz(2.0896623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6733072) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(2.9735273) q[2];
rz(2.7189861) q[3];
sx q[3];
rz(-2.1674619) q[3];
sx q[3];
rz(-0.15338038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733658) q[0];
sx q[0];
rz(-2.234937) q[0];
sx q[0];
rz(-1.9510829) q[0];
rz(-0.52781421) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-0.0091008069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86316696) q[0];
sx q[0];
rz(-0.7896119) q[0];
sx q[0];
rz(-2.4900683) q[0];
rz(-0.78872719) q[2];
sx q[2];
rz(-0.085914748) q[2];
sx q[2];
rz(-0.14363657) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8619662) q[1];
sx q[1];
rz(-1.1512966) q[1];
sx q[1];
rz(2.2395833) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1808257) q[3];
sx q[3];
rz(-2.4159263) q[3];
sx q[3];
rz(-2.8270367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8268383) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(-0.66519386) q[2];
rz(-2.1151755) q[3];
sx q[3];
rz(-2.0668991) q[3];
sx q[3];
rz(2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219532) q[0];
sx q[0];
rz(-2.6060947) q[0];
sx q[0];
rz(2.7435379) q[0];
rz(1.4705426) q[1];
sx q[1];
rz(-2.2649951) q[1];
sx q[1];
rz(-2.3614531) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7791634) q[0];
sx q[0];
rz(-2.0423959) q[0];
sx q[0];
rz(0.40059025) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3684811) q[2];
sx q[2];
rz(-1.8087808) q[2];
sx q[2];
rz(-1.9461103) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4698339) q[1];
sx q[1];
rz(-0.90171725) q[1];
sx q[1];
rz(-1.9356739) q[1];
rz(-pi) q[2];
rz(-2.2557115) q[3];
sx q[3];
rz(-1.8198593) q[3];
sx q[3];
rz(-1.1108171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(-2.9798689) q[2];
rz(3.0297847) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(0.73981729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9991456) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(-2.4430742) q[0];
rz(-1.0922208) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(-3.0879367) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.221401) q[0];
sx q[0];
rz(-2.2760253) q[0];
sx q[0];
rz(2.7979275) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2767792) q[2];
sx q[2];
rz(-0.55241441) q[2];
sx q[2];
rz(2.4557026) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9249869) q[1];
sx q[1];
rz(-2.0408071) q[1];
sx q[1];
rz(1.5108677) q[1];
x q[2];
rz(2.9916006) q[3];
sx q[3];
rz(-2.2120471) q[3];
sx q[3];
rz(-2.7771726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0003537) q[2];
sx q[2];
rz(-2.443479) q[2];
sx q[2];
rz(0.21256438) q[2];
rz(-0.32771787) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(-1.3535961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.245529) q[0];
sx q[0];
rz(-2.0241757) q[0];
sx q[0];
rz(2.8118706) q[0];
rz(-2.3700736) q[1];
sx q[1];
rz(-2.3604269) q[1];
sx q[1];
rz(2.3158997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3701009) q[0];
sx q[0];
rz(-2.2326062) q[0];
sx q[0];
rz(-1.460381) q[0];
rz(-3.0631376) q[2];
sx q[2];
rz(-0.9758853) q[2];
sx q[2];
rz(1.9878146) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4976501) q[1];
sx q[1];
rz(-1.4621328) q[1];
sx q[1];
rz(2.3791299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4903622) q[3];
sx q[3];
rz(-0.31933258) q[3];
sx q[3];
rz(-2.8791219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19994911) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.6597718) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.2991644) q[3];
sx q[3];
rz(0.53946462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.900443) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(-2.7442617) q[0];
rz(-0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(2.1053402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66175211) q[0];
sx q[0];
rz(-1.9755198) q[0];
sx q[0];
rz(-1.5581801) q[0];
rz(-pi) q[1];
rz(-2.7888227) q[2];
sx q[2];
rz(-1.4401541) q[2];
sx q[2];
rz(0.58707159) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4084928) q[1];
sx q[1];
rz(-1.4816135) q[1];
sx q[1];
rz(0.21987889) q[1];
x q[2];
rz(0.37452368) q[3];
sx q[3];
rz(-1.0820884) q[3];
sx q[3];
rz(-0.54218369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.022126023) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(1.3163346) q[2];
rz(-2.8890166) q[3];
sx q[3];
rz(-1.3467237) q[3];
sx q[3];
rz(0.36309567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0040358) q[0];
sx q[0];
rz(-2.3513849) q[0];
sx q[0];
rz(-0.23605119) q[0];
rz(3.0421742) q[1];
sx q[1];
rz(-0.75629083) q[1];
sx q[1];
rz(-1.258446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5262622) q[0];
sx q[0];
rz(-2.1293422) q[0];
sx q[0];
rz(1.6625151) q[0];
rz(-pi) q[1];
x q[1];
rz(2.070932) q[2];
sx q[2];
rz(-0.96958435) q[2];
sx q[2];
rz(-0.76155969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8736781) q[1];
sx q[1];
rz(-2.2482052) q[1];
sx q[1];
rz(-1.4235086) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7332337) q[3];
sx q[3];
rz(-1.1757188) q[3];
sx q[3];
rz(-1.4407106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9749757) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(2.100259) q[2];
rz(-2.5552022) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(-2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36622421) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(-2.1495023) q[1];
sx q[1];
rz(-1.6012259) q[1];
sx q[1];
rz(1.59902) q[1];
rz(-1.3401399) q[2];
sx q[2];
rz(-1.8808095) q[2];
sx q[2];
rz(-1.9145489) q[2];
rz(-0.6057779) q[3];
sx q[3];
rz(-2.45469) q[3];
sx q[3];
rz(-1.5495095) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
