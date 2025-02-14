OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.00081113022) q[0];
sx q[0];
rz(-0.82511628) q[0];
sx q[0];
rz(0.42186475) q[0];
rz(2.3946664) q[1];
sx q[1];
rz(-0.59023017) q[1];
sx q[1];
rz(0.8210558) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79347389) q[0];
sx q[0];
rz(-0.095190053) q[0];
sx q[0];
rz(0.4918672) q[0];
x q[1];
rz(0.73696359) q[2];
sx q[2];
rz(-2.8407477) q[2];
sx q[2];
rz(0.72754169) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6118879) q[1];
sx q[1];
rz(-2.6571353) q[1];
sx q[1];
rz(2.2147708) q[1];
rz(-1.9152568) q[3];
sx q[3];
rz(-1.4446961) q[3];
sx q[3];
rz(-0.9094224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39516285) q[2];
sx q[2];
rz(-0.91150993) q[2];
sx q[2];
rz(2.5968623) q[2];
rz(-1.4970477) q[3];
sx q[3];
rz(-1.112273) q[3];
sx q[3];
rz(1.8018319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4561975) q[0];
sx q[0];
rz(-2.39769) q[0];
sx q[0];
rz(-2.5681382) q[0];
rz(-0.042512976) q[1];
sx q[1];
rz(-1.3923693) q[1];
sx q[1];
rz(-1.858985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1853582) q[0];
sx q[0];
rz(-1.6969674) q[0];
sx q[0];
rz(-1.5915568) q[0];
x q[1];
rz(-1.8268711) q[2];
sx q[2];
rz(-1.8559858) q[2];
sx q[2];
rz(1.8023173) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1309218) q[1];
sx q[1];
rz(-2.0951443) q[1];
sx q[1];
rz(1.4078543) q[1];
rz(-pi) q[2];
rz(2.5124328) q[3];
sx q[3];
rz(-2.0807618) q[3];
sx q[3];
rz(2.8327018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2635228) q[2];
sx q[2];
rz(-1.9922549) q[2];
sx q[2];
rz(2.6980706) q[2];
rz(1.5524607) q[3];
sx q[3];
rz(-0.3905206) q[3];
sx q[3];
rz(0.21903567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1294915) q[0];
sx q[0];
rz(-2.2624367) q[0];
sx q[0];
rz(2.7398859) q[0];
rz(-1.4178287) q[1];
sx q[1];
rz(-2.6938853) q[1];
sx q[1];
rz(2.0254501) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4075466) q[0];
sx q[0];
rz(-2.4067558) q[0];
sx q[0];
rz(-2.6786792) q[0];
x q[1];
rz(-0.87900092) q[2];
sx q[2];
rz(-1.4520116) q[2];
sx q[2];
rz(0.15627651) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.96505839) q[1];
sx q[1];
rz(-1.985834) q[1];
sx q[1];
rz(-1.0227636) q[1];
x q[2];
rz(1.4947555) q[3];
sx q[3];
rz(-2.1180987) q[3];
sx q[3];
rz(0.28441498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22905722) q[2];
sx q[2];
rz(-0.29650649) q[2];
sx q[2];
rz(1.0079481) q[2];
rz(-1.5063162) q[3];
sx q[3];
rz(-1.9624036) q[3];
sx q[3];
rz(2.262871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2849543) q[0];
sx q[0];
rz(-3.1143739) q[0];
sx q[0];
rz(0.83754367) q[0];
rz(-1.8800927) q[1];
sx q[1];
rz(-1.7439758) q[1];
sx q[1];
rz(-1.2417485) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5710302) q[0];
sx q[0];
rz(-1.116291) q[0];
sx q[0];
rz(1.8733371) q[0];
x q[1];
rz(-0.21236944) q[2];
sx q[2];
rz(-2.1687733) q[2];
sx q[2];
rz(-0.47456196) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82363331) q[1];
sx q[1];
rz(-1.1640062) q[1];
sx q[1];
rz(2.645036) q[1];
x q[2];
rz(-2.7697428) q[3];
sx q[3];
rz(-2.1193373) q[3];
sx q[3];
rz(1.8601102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1602829) q[2];
sx q[2];
rz(-0.95429388) q[2];
sx q[2];
rz(2.1295638) q[2];
rz(-2.1435598) q[3];
sx q[3];
rz(-0.73603743) q[3];
sx q[3];
rz(1.6691104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07974682) q[0];
sx q[0];
rz(-2.340402) q[0];
sx q[0];
rz(2.8629942) q[0];
rz(-2.8246236) q[1];
sx q[1];
rz(-1.7520889) q[1];
sx q[1];
rz(-2.680079) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3779031) q[0];
sx q[0];
rz(-0.17556854) q[0];
sx q[0];
rz(-1.2922835) q[0];
x q[1];
rz(-2.1198844) q[2];
sx q[2];
rz(-0.69729348) q[2];
sx q[2];
rz(2.1182228) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0089142) q[1];
sx q[1];
rz(-1.8092833) q[1];
sx q[1];
rz(3.0452252) q[1];
rz(-pi) q[2];
x q[2];
rz(0.096624537) q[3];
sx q[3];
rz(-2.1139675) q[3];
sx q[3];
rz(-2.4718049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3658112) q[2];
sx q[2];
rz(-1.8041939) q[2];
sx q[2];
rz(-0.27654761) q[2];
rz(0.60025275) q[3];
sx q[3];
rz(-2.4252031) q[3];
sx q[3];
rz(-2.6071809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9944331) q[0];
sx q[0];
rz(-2.5614547) q[0];
sx q[0];
rz(-0.68122) q[0];
rz(-2.6874806) q[1];
sx q[1];
rz(-1.157016) q[1];
sx q[1];
rz(-0.62201321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9227064) q[0];
sx q[0];
rz(-1.4758631) q[0];
sx q[0];
rz(-2.9430998) q[0];
x q[1];
rz(-0.63779442) q[2];
sx q[2];
rz(-2.0908815) q[2];
sx q[2];
rz(0.74842036) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4281763) q[1];
sx q[1];
rz(-2.312156) q[1];
sx q[1];
rz(-1.0053289) q[1];
rz(-1.9196904) q[3];
sx q[3];
rz(-1.318228) q[3];
sx q[3];
rz(-0.080772922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5155718) q[2];
sx q[2];
rz(-0.61352789) q[2];
sx q[2];
rz(-0.75508368) q[2];
rz(-0.10609047) q[3];
sx q[3];
rz(-1.2664436) q[3];
sx q[3];
rz(-2.6809926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047091529) q[0];
sx q[0];
rz(-0.61151183) q[0];
sx q[0];
rz(0.066135429) q[0];
rz(-3.0905837) q[1];
sx q[1];
rz(-0.74791932) q[1];
sx q[1];
rz(1.2640818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8569778) q[0];
sx q[0];
rz(-0.98941411) q[0];
sx q[0];
rz(2.7063497) q[0];
rz(-pi) q[1];
rz(-2.4263229) q[2];
sx q[2];
rz(-1.3330622) q[2];
sx q[2];
rz(1.0843074) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5242099) q[1];
sx q[1];
rz(-0.89069378) q[1];
sx q[1];
rz(1.0587803) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43352385) q[3];
sx q[3];
rz(-1.2085087) q[3];
sx q[3];
rz(-1.1690681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2833726) q[2];
sx q[2];
rz(-2.0141352) q[2];
sx q[2];
rz(2.9952725) q[2];
rz(-2.537651) q[3];
sx q[3];
rz(-0.54655176) q[3];
sx q[3];
rz(1.7291791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.716575) q[0];
sx q[0];
rz(-1.1973493) q[0];
sx q[0];
rz(-3.0498411) q[0];
rz(0.07235202) q[1];
sx q[1];
rz(-1.8232583) q[1];
sx q[1];
rz(0.66973698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2429072) q[0];
sx q[0];
rz(-1.7062441) q[0];
sx q[0];
rz(-0.29512914) q[0];
rz(-pi) q[1];
rz(1.4402499) q[2];
sx q[2];
rz(-1.9770844) q[2];
sx q[2];
rz(-0.31995904) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79553855) q[1];
sx q[1];
rz(-0.33949741) q[1];
sx q[1];
rz(-0.44793753) q[1];
x q[2];
rz(-1.1283232) q[3];
sx q[3];
rz(-1.7024346) q[3];
sx q[3];
rz(1.7385755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.79121315) q[2];
sx q[2];
rz(-0.71724856) q[2];
sx q[2];
rz(2.4199602) q[2];
rz(-0.77389884) q[3];
sx q[3];
rz(-0.98267233) q[3];
sx q[3];
rz(-0.53988808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72188193) q[0];
sx q[0];
rz(-1.1996491) q[0];
sx q[0];
rz(0.53959674) q[0];
rz(-0.78060141) q[1];
sx q[1];
rz(-1.2465979) q[1];
sx q[1];
rz(-2.7194729) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2950384) q[0];
sx q[0];
rz(-2.1384412) q[0];
sx q[0];
rz(2.4722383) q[0];
rz(-1.5109748) q[2];
sx q[2];
rz(-1.079353) q[2];
sx q[2];
rz(-2.7697542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.077476689) q[1];
sx q[1];
rz(-1.5327274) q[1];
sx q[1];
rz(-2.0504954) q[1];
x q[2];
rz(-1.6160746) q[3];
sx q[3];
rz(-0.98261925) q[3];
sx q[3];
rz(-2.8578293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.098794) q[2];
sx q[2];
rz(-1.021794) q[2];
sx q[2];
rz(-2.5928024) q[2];
rz(0.15268606) q[3];
sx q[3];
rz(-2.1137674) q[3];
sx q[3];
rz(-2.4299183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3720836) q[0];
sx q[0];
rz(-0.39396572) q[0];
sx q[0];
rz(0.60419303) q[0];
rz(-1.4671885) q[1];
sx q[1];
rz(-1.7908275) q[1];
sx q[1];
rz(-1.9942572) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7580099) q[0];
sx q[0];
rz(-2.1910163) q[0];
sx q[0];
rz(1.0538641) q[0];
x q[1];
rz(0.34530039) q[2];
sx q[2];
rz(-0.7658813) q[2];
sx q[2];
rz(1.6102546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1318329) q[1];
sx q[1];
rz(-1.5554252) q[1];
sx q[1];
rz(-2.6587376) q[1];
x q[2];
rz(1.1710406) q[3];
sx q[3];
rz(-1.2467524) q[3];
sx q[3];
rz(-2.7017665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6818105) q[2];
sx q[2];
rz(-1.1638389) q[2];
sx q[2];
rz(-0.60010827) q[2];
rz(-2.4188304) q[3];
sx q[3];
rz(-0.99812752) q[3];
sx q[3];
rz(-2.7616937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73846524) q[0];
sx q[0];
rz(-1.5140139) q[0];
sx q[0];
rz(-1.4366666) q[0];
rz(-1.5557095) q[1];
sx q[1];
rz(-1.3974421) q[1];
sx q[1];
rz(-0.93346649) q[1];
rz(2.2384833) q[2];
sx q[2];
rz(-1.1739068) q[2];
sx q[2];
rz(0.046233346) q[2];
rz(-1.3550351) q[3];
sx q[3];
rz(-1.329862) q[3];
sx q[3];
rz(-0.71888741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
