OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(-0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(2.1138432) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337024) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(-0.062113751) q[0];
rz(-pi) q[1];
rz(-0.22123863) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(-1.1269119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4029018) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(-2.4163567) q[1];
rz(-pi) q[2];
rz(-0.62698934) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-2.091308) q[2];
rz(1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-2.8438399) q[0];
rz(2.521926) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(2.0334977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28716921) q[0];
sx q[0];
rz(-0.65432917) q[0];
sx q[0];
rz(-1.2229344) q[0];
rz(-pi) q[1];
rz(-2.7116398) q[2];
sx q[2];
rz(-2.1351095) q[2];
sx q[2];
rz(-1.0831837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.83924343) q[1];
sx q[1];
rz(-1.9560555) q[1];
sx q[1];
rz(1.1883931) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3012582) q[3];
sx q[3];
rz(-0.88629913) q[3];
sx q[3];
rz(-1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46253282) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(2.5986824) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(0.96484819) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438) q[0];
sx q[0];
rz(-1.291073) q[0];
sx q[0];
rz(2.8066737) q[0];
rz(-1.5281048) q[2];
sx q[2];
rz(-1.2741158) q[2];
sx q[2];
rz(-1.7745078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0286897) q[1];
sx q[1];
rz(-1.6351055) q[1];
sx q[1];
rz(1.4494004) q[1];
x q[2];
rz(-1.996702) q[3];
sx q[3];
rz(-2.689038) q[3];
sx q[3];
rz(-1.9354613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(-2.0084521) q[2];
rz(-2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(-2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(-1.3765155) q[0];
rz(-2.6230985) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(2.8994697) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74894416) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(1.7783661) q[0];
x q[1];
rz(3.1130586) q[2];
sx q[2];
rz(-2.7118073) q[2];
sx q[2];
rz(2.1536749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.081269216) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(-0.37270765) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4297156) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(-0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(-1.799396) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-0.43911394) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(2.7807996) q[0];
rz(-1.3882673) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(2.0070019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9903588) q[0];
sx q[0];
rz(-0.84999527) q[0];
sx q[0];
rz(-0.6806586) q[0];
x q[1];
rz(-1.1410494) q[2];
sx q[2];
rz(-1.2417214) q[2];
sx q[2];
rz(2.3171901) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6280917) q[1];
sx q[1];
rz(-0.47422945) q[1];
sx q[1];
rz(2.7952816) q[1];
x q[2];
rz(0.43679045) q[3];
sx q[3];
rz(-1.646864) q[3];
sx q[3];
rz(1.2341183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144433) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.1522326) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845247) q[0];
sx q[0];
rz(-3.0355434) q[0];
sx q[0];
rz(-1.4507136) q[0];
x q[1];
rz(2.639159) q[2];
sx q[2];
rz(-1.2049434) q[2];
sx q[2];
rz(-2.1937214) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8730948) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(-0.48392673) q[1];
rz(-pi) q[2];
rz(0.96529393) q[3];
sx q[3];
rz(-2.1552342) q[3];
sx q[3];
rz(2.8402929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(3.0793072) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(2.7468162) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61790723) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(0.54752366) q[0];
x q[1];
rz(2.765352) q[2];
sx q[2];
rz(-1.2708086) q[2];
sx q[2];
rz(-2.8397727) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.863184) q[1];
sx q[1];
rz(-1.8179968) q[1];
sx q[1];
rz(0.43988887) q[1];
rz(-pi) q[2];
rz(1.5338321) q[3];
sx q[3];
rz(-1.7241524) q[3];
sx q[3];
rz(-1.5330029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(2.8685692) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(-2.8295512) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39772314) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.8008908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(-2.7892053) q[0];
rz(-pi) q[1];
rz(1.6216535) q[2];
sx q[2];
rz(-2.2659781) q[2];
sx q[2];
rz(0.67957544) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33412877) q[1];
sx q[1];
rz(-2.6263116) q[1];
sx q[1];
rz(2.3433102) q[1];
x q[2];
rz(-0.31524661) q[3];
sx q[3];
rz(-1.9936221) q[3];
sx q[3];
rz(1.8002312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(2.1179874) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0969365) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(-0.83713371) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9595813) q[0];
sx q[0];
rz(-2.2416229) q[0];
sx q[0];
rz(-1.3120033) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43948549) q[2];
sx q[2];
rz(-2.5186933) q[2];
sx q[2];
rz(2.5650131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3771462) q[1];
sx q[1];
rz(-2.9566777) q[1];
sx q[1];
rz(2.1746693) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0650474) q[3];
sx q[3];
rz(-2.5878083) q[3];
sx q[3];
rz(-1.5050642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(2.0843704) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(-2.1616139) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5532593) q[0];
sx q[0];
rz(-0.13561121) q[0];
sx q[0];
rz(1.9562264) q[0];
rz(-pi) q[1];
rz(2.0699632) q[2];
sx q[2];
rz(-2.0344779) q[2];
sx q[2];
rz(-2.0641363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.92614782) q[1];
sx q[1];
rz(-1.1415592) q[1];
sx q[1];
rz(0.26305671) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2782709) q[3];
sx q[3];
rz(-2.8046126) q[3];
sx q[3];
rz(-1.8388336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86823157) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(2.3095619) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(1.8503415) q[2];
sx q[2];
rz(-0.57689473) q[2];
sx q[2];
rz(2.9070791) q[2];
rz(2.5084393) q[3];
sx q[3];
rz(-2.5452151) q[3];
sx q[3];
rz(2.7635318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
