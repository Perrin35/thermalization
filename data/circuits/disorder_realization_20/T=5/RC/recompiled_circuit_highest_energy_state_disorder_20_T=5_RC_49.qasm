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
rz(0.19658495) q[0];
sx q[0];
rz(-2.6994446) q[0];
sx q[0];
rz(-0.86139876) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(-1.8169401) q[1];
sx q[1];
rz(3.1060001) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4342049) q[0];
sx q[0];
rz(-1.905283) q[0];
sx q[0];
rz(-1.3637961) q[0];
rz(2.8096417) q[2];
sx q[2];
rz(-2.8185726) q[2];
sx q[2];
rz(-2.2145751) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.030499996) q[1];
sx q[1];
rz(-0.72734088) q[1];
sx q[1];
rz(-0.05383454) q[1];
rz(-3.1236105) q[3];
sx q[3];
rz(-0.89089046) q[3];
sx q[3];
rz(2.3462092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5888136) q[2];
sx q[2];
rz(-1.3014883) q[2];
sx q[2];
rz(-1.0944875) q[2];
rz(-1.6257446) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(2.998013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0482408) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(-0.37407237) q[0];
rz(0.85809842) q[1];
sx q[1];
rz(-2.2007807) q[1];
sx q[1];
rz(2.0657952) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5257403) q[0];
sx q[0];
rz(-1.6739556) q[0];
sx q[0];
rz(3.0391058) q[0];
rz(-pi) q[1];
rz(-2.7611742) q[2];
sx q[2];
rz(-0.58092344) q[2];
sx q[2];
rz(-0.38512938) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0829216) q[1];
sx q[1];
rz(-1.6870912) q[1];
sx q[1];
rz(0.40177675) q[1];
rz(0.91173633) q[3];
sx q[3];
rz(-1.0056595) q[3];
sx q[3];
rz(2.679058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51521236) q[2];
sx q[2];
rz(-2.7084646) q[2];
sx q[2];
rz(-2.058378) q[2];
rz(1.8488688) q[3];
sx q[3];
rz(-2.0079398) q[3];
sx q[3];
rz(-2.2129272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2556297) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-2.6166925) q[0];
rz(1.6820172) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(-0.14579138) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4065259) q[0];
sx q[0];
rz(-0.13938381) q[0];
sx q[0];
rz(2.4295619) q[0];
rz(-2.2810788) q[2];
sx q[2];
rz(-0.73627055) q[2];
sx q[2];
rz(-2.3585573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5151095) q[1];
sx q[1];
rz(-1.4806387) q[1];
sx q[1];
rz(0.43403352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94696857) q[3];
sx q[3];
rz(-0.83443975) q[3];
sx q[3];
rz(2.8744227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8968481) q[2];
sx q[2];
rz(-1.9795828) q[2];
sx q[2];
rz(-2.9761918) q[2];
rz(-2.2798955) q[3];
sx q[3];
rz(-2.4498037) q[3];
sx q[3];
rz(0.19680463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621314) q[0];
sx q[0];
rz(-0.52424163) q[0];
sx q[0];
rz(-2.4133546) q[0];
rz(-2.8184452) q[1];
sx q[1];
rz(-2.1073982) q[1];
sx q[1];
rz(2.1023777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6704039) q[0];
sx q[0];
rz(-1.146933) q[0];
sx q[0];
rz(-0.43760646) q[0];
rz(-1.9437691) q[2];
sx q[2];
rz(-0.67016888) q[2];
sx q[2];
rz(-3.016346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2715214) q[1];
sx q[1];
rz(-2.0539375) q[1];
sx q[1];
rz(1.6752233) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3023241) q[3];
sx q[3];
rz(-0.89863649) q[3];
sx q[3];
rz(1.9306926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49179658) q[2];
sx q[2];
rz(-1.0504664) q[2];
sx q[2];
rz(0.51898471) q[2];
rz(-1.0745878) q[3];
sx q[3];
rz(-1.855775) q[3];
sx q[3];
rz(0.49526596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.062926) q[0];
sx q[0];
rz(-0.92200297) q[0];
sx q[0];
rz(2.9700188) q[0];
rz(1.0942787) q[1];
sx q[1];
rz(-1.8405874) q[1];
sx q[1];
rz(2.9069854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6162524) q[0];
sx q[0];
rz(-1.8484637) q[0];
sx q[0];
rz(0.13292952) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.094298074) q[2];
sx q[2];
rz(-2.6324632) q[2];
sx q[2];
rz(-1.5056339) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.460894) q[1];
sx q[1];
rz(-2.1354224) q[1];
sx q[1];
rz(2.2959397) q[1];
x q[2];
rz(-2.4314546) q[3];
sx q[3];
rz(-1.8605068) q[3];
sx q[3];
rz(-2.3848145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2285063) q[2];
sx q[2];
rz(-1.191491) q[2];
sx q[2];
rz(-1.7388434) q[2];
rz(2.5168822) q[3];
sx q[3];
rz(-0.96840817) q[3];
sx q[3];
rz(-1.3431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6828571) q[0];
sx q[0];
rz(-3.1338437) q[0];
sx q[0];
rz(0.32753456) q[0];
rz(0.028118357) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(-2.1498674) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018143749) q[0];
sx q[0];
rz(-0.94052343) q[0];
sx q[0];
rz(-0.92395497) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4621134) q[2];
sx q[2];
rz(-1.7955488) q[2];
sx q[2];
rz(1.0959742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8918482) q[1];
sx q[1];
rz(-1.5288903) q[1];
sx q[1];
rz(-0.81467198) q[1];
rz(-pi) q[2];
x q[2];
rz(3.087849) q[3];
sx q[3];
rz(-2.6702767) q[3];
sx q[3];
rz(2.0278679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.907054) q[2];
sx q[2];
rz(-1.7639284) q[2];
sx q[2];
rz(-1.6592337) q[2];
rz(1.0659069) q[3];
sx q[3];
rz(-1.1803455) q[3];
sx q[3];
rz(-1.0970241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3858353) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(0.29543153) q[0];
rz(2.9040728) q[1];
sx q[1];
rz(-0.93370456) q[1];
sx q[1];
rz(-2.3942153) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.462446) q[0];
sx q[0];
rz(-2.0751795) q[0];
sx q[0];
rz(-2.4997968) q[0];
rz(-2.4557023) q[2];
sx q[2];
rz(-0.7524366) q[2];
sx q[2];
rz(2.8519832) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.805757) q[1];
sx q[1];
rz(-1.079396) q[1];
sx q[1];
rz(-0.42224981) q[1];
rz(-pi) q[2];
rz(2.564424) q[3];
sx q[3];
rz(-0.62466972) q[3];
sx q[3];
rz(0.28055252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4733009) q[2];
sx q[2];
rz(-1.95582) q[2];
sx q[2];
rz(0.2549003) q[2];
rz(-2.9292246) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(-3.1406241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85106987) q[0];
sx q[0];
rz(-1.2238598) q[0];
sx q[0];
rz(3.0176924) q[0];
rz(-0.87521416) q[1];
sx q[1];
rz(-0.83299914) q[1];
sx q[1];
rz(0.55878729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13128114) q[0];
sx q[0];
rz(-1.049066) q[0];
sx q[0];
rz(-2.8521796) q[0];
x q[1];
rz(-2.0622425) q[2];
sx q[2];
rz(-0.55365593) q[2];
sx q[2];
rz(1.4244256) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3271823) q[1];
sx q[1];
rz(-2.5375527) q[1];
sx q[1];
rz(-2.6534897) q[1];
x q[2];
rz(0.41591538) q[3];
sx q[3];
rz(-0.97813604) q[3];
sx q[3];
rz(2.376568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.045804068) q[2];
sx q[2];
rz(-0.49171058) q[2];
sx q[2];
rz(0.72074214) q[2];
rz(0.36302429) q[3];
sx q[3];
rz(-1.9178773) q[3];
sx q[3];
rz(1.629841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95661288) q[0];
sx q[0];
rz(-0.19421254) q[0];
sx q[0];
rz(-0.089056253) q[0];
rz(-0.36674276) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(1.4595703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3066481) q[0];
sx q[0];
rz(-0.95150009) q[0];
sx q[0];
rz(2.8602116) q[0];
x q[1];
rz(-2.4581681) q[2];
sx q[2];
rz(-2.0439842) q[2];
sx q[2];
rz(0.13147182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1398581) q[1];
sx q[1];
rz(-2.325104) q[1];
sx q[1];
rz(-2.0862608) q[1];
rz(-pi) q[2];
rz(-0.59883786) q[3];
sx q[3];
rz(-1.2019079) q[3];
sx q[3];
rz(1.9112596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2064994) q[2];
sx q[2];
rz(-1.3672028) q[2];
sx q[2];
rz(1.2500786) q[2];
rz(-1.2262723) q[3];
sx q[3];
rz(-0.63392249) q[3];
sx q[3];
rz(-2.5026076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298237) q[0];
sx q[0];
rz(-2.0084232) q[0];
sx q[0];
rz(2.0015707) q[0];
rz(0.58311588) q[1];
sx q[1];
rz(-1.6551599) q[1];
sx q[1];
rz(-0.66782943) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.697628) q[0];
sx q[0];
rz(-0.85845843) q[0];
sx q[0];
rz(-1.1348073) q[0];
x q[1];
rz(1.2905707) q[2];
sx q[2];
rz(-2.4200768) q[2];
sx q[2];
rz(-2.1211598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7857504) q[1];
sx q[1];
rz(-1.3763345) q[1];
sx q[1];
rz(2.0114312) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79270122) q[3];
sx q[3];
rz(-1.2374733) q[3];
sx q[3];
rz(1.147721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84393152) q[2];
sx q[2];
rz(-0.49489489) q[2];
sx q[2];
rz(-2.3893791) q[2];
rz(-0.080605896) q[3];
sx q[3];
rz(-1.7919431) q[3];
sx q[3];
rz(-2.5940671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3362296) q[0];
sx q[0];
rz(-1.5267876) q[0];
sx q[0];
rz(3.0841854) q[0];
rz(-0.84857955) q[1];
sx q[1];
rz(-2.4129557) q[1];
sx q[1];
rz(1.6853263) q[1];
rz(0.91114974) q[2];
sx q[2];
rz(-1.0187889) q[2];
sx q[2];
rz(2.7762882) q[2];
rz(2.6808394) q[3];
sx q[3];
rz(-1.2450153) q[3];
sx q[3];
rz(1.0871441) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
