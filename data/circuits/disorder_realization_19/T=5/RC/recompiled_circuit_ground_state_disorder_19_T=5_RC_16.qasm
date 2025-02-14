OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8813397) q[0];
sx q[0];
rz(-0.94085675) q[0];
sx q[0];
rz(2.9139304) q[0];
rz(0.28336278) q[1];
sx q[1];
rz(3.5609666) q[1];
sx q[1];
rz(11.279439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.697098) q[0];
sx q[0];
rz(-0.47045204) q[0];
sx q[0];
rz(1.4509499) q[0];
rz(-0.46704328) q[2];
sx q[2];
rz(-1.1640295) q[2];
sx q[2];
rz(-1.9860991) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4408001) q[1];
sx q[1];
rz(-1.9682574) q[1];
sx q[1];
rz(2.1093919) q[1];
x q[2];
rz(-0.6391949) q[3];
sx q[3];
rz(-1.3167644) q[3];
sx q[3];
rz(-0.81786903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2490354) q[2];
sx q[2];
rz(-1.0113357) q[2];
sx q[2];
rz(-0.91157836) q[2];
rz(-2.3859731) q[3];
sx q[3];
rz(-2.8204462) q[3];
sx q[3];
rz(0.49629456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3590473) q[0];
sx q[0];
rz(-1.3107212) q[0];
sx q[0];
rz(-1.1388592) q[0];
rz(0.62659872) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(2.0638594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12909938) q[0];
sx q[0];
rz(-1.9097304) q[0];
sx q[0];
rz(0.61130543) q[0];
rz(-0.90684129) q[2];
sx q[2];
rz(-2.5996947) q[2];
sx q[2];
rz(-0.52846891) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6328395) q[1];
sx q[1];
rz(-1.6378073) q[1];
sx q[1];
rz(0.26439039) q[1];
rz(2.7116345) q[3];
sx q[3];
rz(-0.89498496) q[3];
sx q[3];
rz(1.2682008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44800147) q[2];
sx q[2];
rz(-2.3841264) q[2];
sx q[2];
rz(-0.21270154) q[2];
rz(-1.5851783) q[3];
sx q[3];
rz(-2.3978265) q[3];
sx q[3];
rz(2.0785544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146516) q[0];
sx q[0];
rz(-0.26697049) q[0];
sx q[0];
rz(-0.84685999) q[0];
rz(0.2521387) q[1];
sx q[1];
rz(-0.82770258) q[1];
sx q[1];
rz(-0.65021461) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51948773) q[0];
sx q[0];
rz(-1.5267924) q[0];
sx q[0];
rz(-0.33885886) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6619307) q[2];
sx q[2];
rz(-2.227894) q[2];
sx q[2];
rz(-1.3160694) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.099459186) q[1];
sx q[1];
rz(-0.317527) q[1];
sx q[1];
rz(1.2088163) q[1];
rz(0.80649431) q[3];
sx q[3];
rz(-1.8915081) q[3];
sx q[3];
rz(2.210828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.83027679) q[2];
sx q[2];
rz(-0.69861424) q[2];
sx q[2];
rz(2.1777731) q[2];
rz(1.7802995) q[3];
sx q[3];
rz(-2.6908974) q[3];
sx q[3];
rz(3.0986339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394102) q[0];
sx q[0];
rz(-2.9189411) q[0];
sx q[0];
rz(-1.0182925) q[0];
rz(-2.2564383) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(0.28908602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0196298) q[0];
sx q[0];
rz(-1.2696854) q[0];
sx q[0];
rz(0.34261468) q[0];
rz(1.1491597) q[2];
sx q[2];
rz(-2.1648063) q[2];
sx q[2];
rz(-0.27911148) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.30384597) q[1];
sx q[1];
rz(-1.0582542) q[1];
sx q[1];
rz(2.9339177) q[1];
x q[2];
rz(2.9226275) q[3];
sx q[3];
rz(-1.4462399) q[3];
sx q[3];
rz(-1.7809465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3503512) q[2];
sx q[2];
rz(-1.505932) q[2];
sx q[2];
rz(-1.852847) q[2];
rz(-2.5617803) q[3];
sx q[3];
rz(-2.6668187) q[3];
sx q[3];
rz(-0.60540664) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30699214) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(-2.8173764) q[0];
rz(3.1027555) q[1];
sx q[1];
rz(-2.4034998) q[1];
sx q[1];
rz(-2.0447581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997579) q[0];
sx q[0];
rz(-1.2909162) q[0];
sx q[0];
rz(1.8246234) q[0];
x q[1];
rz(-2.1514405) q[2];
sx q[2];
rz(-1.6564995) q[2];
sx q[2];
rz(3.1243968) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.005744) q[1];
sx q[1];
rz(-1.6127819) q[1];
sx q[1];
rz(-1.2688925) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9911258) q[3];
sx q[3];
rz(-1.655135) q[3];
sx q[3];
rz(-1.4140417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0754806) q[2];
sx q[2];
rz(-2.8754063) q[2];
sx q[2];
rz(-0.038662635) q[2];
rz(-1.4085116) q[3];
sx q[3];
rz(-1.6950636) q[3];
sx q[3];
rz(-2.8601638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51327813) q[0];
sx q[0];
rz(-2.3230041) q[0];
sx q[0];
rz(1.0191089) q[0];
rz(-2.6844773) q[1];
sx q[1];
rz(-2.8725084) q[1];
sx q[1];
rz(1.7105182) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015648677) q[0];
sx q[0];
rz(-1.5971184) q[0];
sx q[0];
rz(-0.20080815) q[0];
x q[1];
rz(0.25361714) q[2];
sx q[2];
rz(-2.1448958) q[2];
sx q[2];
rz(0.92437896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3574111) q[1];
sx q[1];
rz(-0.37115449) q[1];
sx q[1];
rz(1.3264912) q[1];
rz(-pi) q[2];
rz(-2.2912628) q[3];
sx q[3];
rz(-1.2652694) q[3];
sx q[3];
rz(3.1269493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9223601) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(-0.57979453) q[2];
rz(-2.842105) q[3];
sx q[3];
rz(-0.53286415) q[3];
sx q[3];
rz(-1.2286435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806468) q[0];
sx q[0];
rz(-1.6609284) q[0];
sx q[0];
rz(0.10375599) q[0];
rz(0.62458986) q[1];
sx q[1];
rz(-0.92085212) q[1];
sx q[1];
rz(1.3821028) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64829761) q[0];
sx q[0];
rz(-1.7837423) q[0];
sx q[0];
rz(-2.2981724) q[0];
rz(-pi) q[1];
rz(-0.29989179) q[2];
sx q[2];
rz(-1.4647398) q[2];
sx q[2];
rz(0.40115717) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8727726) q[1];
sx q[1];
rz(-1.8389069) q[1];
sx q[1];
rz(1.0467749) q[1];
rz(-pi) q[2];
rz(0.66094605) q[3];
sx q[3];
rz(-1.3437643) q[3];
sx q[3];
rz(1.1477136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3524807) q[2];
sx q[2];
rz(-1.4752957) q[2];
sx q[2];
rz(2.2769807) q[2];
rz(2.9027446) q[3];
sx q[3];
rz(-2.3819203) q[3];
sx q[3];
rz(2.4281832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9637941) q[0];
sx q[0];
rz(-0.56121427) q[0];
sx q[0];
rz(-2.906565) q[0];
rz(2.4759953) q[1];
sx q[1];
rz(-2.0033629) q[1];
sx q[1];
rz(-1.4733018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99645961) q[0];
sx q[0];
rz(-2.1204824) q[0];
sx q[0];
rz(0.24497801) q[0];
x q[1];
rz(-0.46730475) q[2];
sx q[2];
rz(-1.4479785) q[2];
sx q[2];
rz(-0.87087028) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9884432) q[1];
sx q[1];
rz(-1.4979939) q[1];
sx q[1];
rz(2.9576587) q[1];
x q[2];
rz(3.0489122) q[3];
sx q[3];
rz(-2.6887935) q[3];
sx q[3];
rz(-1.6860675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9792446) q[2];
sx q[2];
rz(-0.20077106) q[2];
sx q[2];
rz(2.6806504) q[2];
rz(-1.0848684) q[3];
sx q[3];
rz(-0.74551398) q[3];
sx q[3];
rz(-2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.064086176) q[0];
sx q[0];
rz(-2.6106847) q[0];
sx q[0];
rz(2.532646) q[0];
rz(2.5755836) q[1];
sx q[1];
rz(-1.1606263) q[1];
sx q[1];
rz(1.2014679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.214819) q[0];
sx q[0];
rz(-1.578971) q[0];
sx q[0];
rz(-1.4352485) q[0];
rz(1.1913774) q[2];
sx q[2];
rz(-2.622329) q[2];
sx q[2];
rz(1.0644703) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4200099) q[1];
sx q[1];
rz(-0.8045361) q[1];
sx q[1];
rz(-1.4415506) q[1];
rz(-0.44916656) q[3];
sx q[3];
rz(-1.2954062) q[3];
sx q[3];
rz(-0.034253293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0386049) q[2];
sx q[2];
rz(-2.0709585) q[2];
sx q[2];
rz(-1.9496244) q[2];
rz(-2.1206756) q[3];
sx q[3];
rz(-2.9396785) q[3];
sx q[3];
rz(-1.5405704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9182619) q[0];
sx q[0];
rz(-0.20063618) q[0];
sx q[0];
rz(-0.9675135) q[0];
rz(-1.7009023) q[1];
sx q[1];
rz(-2.7099425) q[1];
sx q[1];
rz(-0.38223019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0456881) q[0];
sx q[0];
rz(-1.75215) q[0];
sx q[0];
rz(-2.1482094) q[0];
rz(-pi) q[1];
rz(-0.43469825) q[2];
sx q[2];
rz(-1.189173) q[2];
sx q[2];
rz(-1.3023072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2652685) q[1];
sx q[1];
rz(-2.3493715) q[1];
sx q[1];
rz(1.4033615) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1145544) q[3];
sx q[3];
rz(-0.54618764) q[3];
sx q[3];
rz(-2.9423713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.050015673) q[2];
sx q[2];
rz(-0.86805582) q[2];
sx q[2];
rz(-0.78563219) q[2];
rz(-0.026963726) q[3];
sx q[3];
rz(-1.5188768) q[3];
sx q[3];
rz(0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.3383269) q[0];
sx q[0];
rz(-1.6852408) q[0];
sx q[0];
rz(2.2483873) q[0];
rz(0.66435736) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(-2.4128466) q[2];
sx q[2];
rz(-1.3441636) q[2];
sx q[2];
rz(-2.8580479) q[2];
rz(2.5369) q[3];
sx q[3];
rz(-1.5786377) q[3];
sx q[3];
rz(-0.64463401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
