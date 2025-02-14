OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26025298) q[0];
sx q[0];
rz(4.0824494) q[0];
sx q[0];
rz(9.6524402) q[0];
rz(-2.8582299) q[1];
sx q[1];
rz(-0.41937399) q[1];
sx q[1];
rz(-1.8546606) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5628107) q[0];
sx q[0];
rz(-2.037604) q[0];
sx q[0];
rz(-0.060725529) q[0];
rz(0.76331933) q[2];
sx q[2];
rz(-0.60930353) q[2];
sx q[2];
rz(-2.0610025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4408001) q[1];
sx q[1];
rz(-1.1733353) q[1];
sx q[1];
rz(-2.1093919) q[1];
x q[2];
rz(0.41051045) q[3];
sx q[3];
rz(-0.68119739) q[3];
sx q[3];
rz(2.714701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2490354) q[2];
sx q[2];
rz(-1.0113357) q[2];
sx q[2];
rz(2.2300143) q[2];
rz(-0.75561953) q[3];
sx q[3];
rz(-0.32114649) q[3];
sx q[3];
rz(0.49629456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78254533) q[0];
sx q[0];
rz(-1.3107212) q[0];
sx q[0];
rz(1.1388592) q[0];
rz(2.5149939) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(1.0777333) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12909938) q[0];
sx q[0];
rz(-1.2318622) q[0];
sx q[0];
rz(-0.61130543) q[0];
rz(-pi) q[1];
rz(-0.35524345) q[2];
sx q[2];
rz(-1.9890824) q[2];
sx q[2];
rz(0.21165161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0976757) q[1];
sx q[1];
rz(-1.8345791) q[1];
sx q[1];
rz(1.640212) q[1];
x q[2];
rz(2.7116345) q[3];
sx q[3];
rz(-0.89498496) q[3];
sx q[3];
rz(1.2682008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6935912) q[2];
sx q[2];
rz(-0.75746626) q[2];
sx q[2];
rz(2.9288911) q[2];
rz(-1.5851783) q[3];
sx q[3];
rz(-0.74376619) q[3];
sx q[3];
rz(-2.0785544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146516) q[0];
sx q[0];
rz(-2.8746222) q[0];
sx q[0];
rz(-2.2947327) q[0];
rz(0.2521387) q[1];
sx q[1];
rz(-2.3138901) q[1];
sx q[1];
rz(-2.491378) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0747809) q[0];
sx q[0];
rz(-1.2322786) q[0];
sx q[0];
rz(-1.5241429) q[0];
x q[1];
rz(-0.89805897) q[2];
sx q[2];
rz(-2.2453547) q[2];
sx q[2];
rz(-2.7314699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47874988) q[1];
sx q[1];
rz(-1.2744941) q[1];
sx q[1];
rz(0.11586145) q[1];
rz(2.0183205) q[3];
sx q[3];
rz(-2.3254804) q[3];
sx q[3];
rz(2.1838674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83027679) q[2];
sx q[2];
rz(-0.69861424) q[2];
sx q[2];
rz(2.1777731) q[2];
rz(1.3612932) q[3];
sx q[3];
rz(-0.45069525) q[3];
sx q[3];
rz(-0.042958766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90218246) q[0];
sx q[0];
rz(-2.9189411) q[0];
sx q[0];
rz(1.0182925) q[0];
rz(2.2564383) q[1];
sx q[1];
rz(-0.2489018) q[1];
sx q[1];
rz(-2.8525066) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24482306) q[0];
sx q[0];
rz(-0.45216783) q[0];
sx q[0];
rz(-0.74613476) q[0];
rz(-1.1491597) q[2];
sx q[2];
rz(-2.1648063) q[2];
sx q[2];
rz(0.27911148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4319181) q[1];
sx q[1];
rz(-0.54952114) q[1];
sx q[1];
rz(-1.9220244) q[1];
rz(1.4432256) q[3];
sx q[3];
rz(-1.3535548) q[3];
sx q[3];
rz(-2.9590817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79124147) q[2];
sx q[2];
rz(-1.505932) q[2];
sx q[2];
rz(1.852847) q[2];
rz(2.5617803) q[3];
sx q[3];
rz(-0.47477397) q[3];
sx q[3];
rz(-0.60540664) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8346005) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(0.32421625) q[0];
rz(0.038837198) q[1];
sx q[1];
rz(-0.7380929) q[1];
sx q[1];
rz(-2.0447581) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997579) q[0];
sx q[0];
rz(-1.8506764) q[0];
sx q[0];
rz(-1.8246234) q[0];
x q[1];
rz(-2.1514405) q[2];
sx q[2];
rz(-1.6564995) q[2];
sx q[2];
rz(3.1243968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.13584863) q[1];
sx q[1];
rz(-1.5288108) q[1];
sx q[1];
rz(1.2688925) q[1];
rz(-pi) q[2];
rz(-0.092336307) q[3];
sx q[3];
rz(-1.9895377) q[3];
sx q[3];
rz(-0.11912042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0754806) q[2];
sx q[2];
rz(-2.8754063) q[2];
sx q[2];
rz(3.10293) q[2];
rz(1.4085116) q[3];
sx q[3];
rz(-1.446529) q[3];
sx q[3];
rz(0.2814289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51327813) q[0];
sx q[0];
rz(-0.81858855) q[0];
sx q[0];
rz(-2.1224838) q[0];
rz(2.6844773) q[1];
sx q[1];
rz(-2.8725084) q[1];
sx q[1];
rz(-1.7105182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5810878) q[0];
sx q[0];
rz(-1.771534) q[0];
sx q[0];
rz(-1.597658) q[0];
rz(-pi) q[1];
rz(0.25361714) q[2];
sx q[2];
rz(-2.1448958) q[2];
sx q[2];
rz(0.92437896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3574111) q[1];
sx q[1];
rz(-0.37115449) q[1];
sx q[1];
rz(1.3264912) q[1];
rz(-pi) q[2];
rz(2.2912628) q[3];
sx q[3];
rz(-1.2652694) q[3];
sx q[3];
rz(0.014643365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9223601) q[2];
sx q[2];
rz(-1.3693501) q[2];
sx q[2];
rz(2.5617981) q[2];
rz(2.842105) q[3];
sx q[3];
rz(-2.6087285) q[3];
sx q[3];
rz(-1.2286435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806468) q[0];
sx q[0];
rz(-1.6609284) q[0];
sx q[0];
rz(-3.0378367) q[0];
rz(2.5170028) q[1];
sx q[1];
rz(-0.92085212) q[1];
sx q[1];
rz(-1.3821028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860024) q[0];
sx q[0];
rz(-0.75240032) q[0];
sx q[0];
rz(1.2563906) q[0];
rz(0.34587282) q[2];
sx q[2];
rz(-2.8240339) q[2];
sx q[2];
rz(1.6421184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8727726) q[1];
sx q[1];
rz(-1.8389069) q[1];
sx q[1];
rz(2.0948177) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66094605) q[3];
sx q[3];
rz(-1.7978284) q[3];
sx q[3];
rz(-1.1477136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.78911191) q[2];
sx q[2];
rz(-1.4752957) q[2];
sx q[2];
rz(2.2769807) q[2];
rz(0.23884808) q[3];
sx q[3];
rz(-0.75967234) q[3];
sx q[3];
rz(2.4281832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1777986) q[0];
sx q[0];
rz(-0.56121427) q[0];
sx q[0];
rz(-2.906565) q[0];
rz(-0.66559732) q[1];
sx q[1];
rz(-1.1382297) q[1];
sx q[1];
rz(-1.6682909) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6971254) q[0];
sx q[0];
rz(-1.3624862) q[0];
sx q[0];
rz(1.0075157) q[0];
x q[1];
rz(0.46730475) q[2];
sx q[2];
rz(-1.6936142) q[2];
sx q[2];
rz(-0.87087028) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.044933407) q[1];
sx q[1];
rz(-0.19766624) q[1];
sx q[1];
rz(-2.7621619) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0489122) q[3];
sx q[3];
rz(-2.6887935) q[3];
sx q[3];
rz(-1.4555252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9792446) q[2];
sx q[2];
rz(-0.20077106) q[2];
sx q[2];
rz(0.46094224) q[2];
rz(-2.0567242) q[3];
sx q[3];
rz(-2.3960787) q[3];
sx q[3];
rz(-2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-3.0775065) q[0];
sx q[0];
rz(-2.6106847) q[0];
sx q[0];
rz(-2.532646) q[0];
rz(0.56600904) q[1];
sx q[1];
rz(-1.9809664) q[1];
sx q[1];
rz(-1.9401248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5574359) q[0];
sx q[0];
rz(-0.13579255) q[0];
sx q[0];
rz(-1.5103755) q[0];
x q[1];
rz(-1.0827093) q[2];
sx q[2];
rz(-1.3859473) q[2];
sx q[2];
rz(-0.83959296) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9068599) q[1];
sx q[1];
rz(-0.77488778) q[1];
sx q[1];
rz(-3.0084684) q[1];
rz(0.44916656) q[3];
sx q[3];
rz(-1.2954062) q[3];
sx q[3];
rz(0.034253293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0386049) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(1.1919682) q[2];
rz(-1.0209171) q[3];
sx q[3];
rz(-0.20191419) q[3];
sx q[3];
rz(-1.5405704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22333071) q[0];
sx q[0];
rz(-2.9409565) q[0];
sx q[0];
rz(-2.1740792) q[0];
rz(-1.4406904) q[1];
sx q[1];
rz(-2.7099425) q[1];
sx q[1];
rz(0.38223019) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8865693) q[0];
sx q[0];
rz(-0.60211997) q[0];
sx q[0];
rz(1.894879) q[0];
rz(-pi) q[1];
x q[1];
rz(2.380314) q[2];
sx q[2];
rz(-0.57028162) q[2];
sx q[2];
rz(-0.40752652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5014508) q[1];
sx q[1];
rz(-0.79265187) q[1];
sx q[1];
rz(-0.16736302) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1145544) q[3];
sx q[3];
rz(-2.595405) q[3];
sx q[3];
rz(2.9423713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.091577) q[2];
sx q[2];
rz(-2.2735368) q[2];
sx q[2];
rz(2.3559605) q[2];
rz(-3.1146289) q[3];
sx q[3];
rz(-1.6227159) q[3];
sx q[3];
rz(0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80326573) q[0];
sx q[0];
rz(-1.6852408) q[0];
sx q[0];
rz(2.2483873) q[0];
rz(-2.4772353) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(-2.4128466) q[2];
sx q[2];
rz(-1.3441636) q[2];
sx q[2];
rz(-2.8580479) q[2];
rz(0.01379225) q[3];
sx q[3];
rz(-0.60473718) q[3];
sx q[3];
rz(-2.2040839) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
