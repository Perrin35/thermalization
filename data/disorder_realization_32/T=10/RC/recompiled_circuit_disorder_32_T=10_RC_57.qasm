OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(2.0523235) q[0];
sx q[0];
rz(9.2637445) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8344603) q[0];
sx q[0];
rz(-1.4559329) q[0];
sx q[0];
rz(0.97712028) q[0];
rz(-pi) q[1];
rz(2.6167294) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(2.2847069) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1094184) q[1];
sx q[1];
rz(-0.86583455) q[1];
sx q[1];
rz(1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6707889) q[3];
sx q[3];
rz(-1.1649719) q[3];
sx q[3];
rz(-1.0909181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(2.2709133) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(1.5354935) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664415) q[0];
sx q[0];
rz(-1.6834714) q[0];
sx q[0];
rz(-2.9215647) q[0];
rz(-2.8385542) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(2.9532202) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9175375) q[1];
sx q[1];
rz(-1.1174669) q[1];
sx q[1];
rz(-0.097150306) q[1];
rz(2.9868449) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(-0.89150067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(1.0822901) q[2];
rz(1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(-0.043116365) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(-1.9972237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17830081) q[0];
sx q[0];
rz(-0.68094567) q[0];
sx q[0];
rz(-0.91239022) q[0];
rz(-pi) q[1];
rz(-1.794533) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(-0.70659107) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.80013393) q[1];
sx q[1];
rz(-1.1143436) q[1];
sx q[1];
rz(1.9564499) q[1];
rz(-pi) q[2];
rz(-3.1312917) q[3];
sx q[3];
rz(-2.8539538) q[3];
sx q[3];
rz(-1.079263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.959789) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(-0.71587193) q[0];
rz(-2.8158358) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(2.148927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0621322) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(-2.604565) q[0];
rz(-pi) q[1];
rz(-2.2067581) q[2];
sx q[2];
rz(-1.0996498) q[2];
sx q[2];
rz(-0.83540321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.042295481) q[1];
sx q[1];
rz(-0.18810939) q[1];
sx q[1];
rz(0.52109615) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99289258) q[3];
sx q[3];
rz(-1.143647) q[3];
sx q[3];
rz(-0.072673365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(0.76675057) q[0];
rz(-1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-2.7899182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7809217) q[0];
sx q[0];
rz(-0.70510266) q[0];
sx q[0];
rz(-0.4522001) q[0];
x q[1];
rz(-2.211116) q[2];
sx q[2];
rz(-0.74379197) q[2];
sx q[2];
rz(-3.004068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1406447) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(2.7468202) q[1];
rz(-2.8204927) q[3];
sx q[3];
rz(-2.0541111) q[3];
sx q[3];
rz(-1.6881642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(0.53331214) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(-2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(1.0995964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27605155) q[0];
sx q[0];
rz(-1.2035032) q[0];
sx q[0];
rz(-2.0588576) q[0];
rz(-pi) q[1];
rz(2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(-2.1649233) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53193608) q[1];
sx q[1];
rz(-1.3927288) q[1];
sx q[1];
rz(2.460536) q[1];
rz(-pi) q[2];
rz(-1.908329) q[3];
sx q[3];
rz(-0.81157717) q[3];
sx q[3];
rz(2.8525762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(-2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-2.506536) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866078) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(2.9842916) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3646506) q[0];
sx q[0];
rz(-1.6645263) q[0];
sx q[0];
rz(-2.4073497) q[0];
rz(-pi) q[1];
rz(0.21643164) q[2];
sx q[2];
rz(-0.84605233) q[2];
sx q[2];
rz(3.016678) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3992608) q[1];
sx q[1];
rz(-0.6472339) q[1];
sx q[1];
rz(-1.2433692) q[1];
x q[2];
rz(0.62742426) q[3];
sx q[3];
rz(-1.3324696) q[3];
sx q[3];
rz(2.0605007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9592231) q[2];
sx q[2];
rz(-0.59252858) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4956932) q[0];
sx q[0];
rz(-2.2024787) q[0];
sx q[0];
rz(0.010852531) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(-0.29327926) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0115259) q[0];
sx q[0];
rz(-1.5648601) q[0];
sx q[0];
rz(-2.5991711) q[0];
rz(-pi) q[1];
x q[1];
rz(1.79004) q[2];
sx q[2];
rz(-0.65537894) q[2];
sx q[2];
rz(-3.0041681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7640904) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(2.3368895) q[1];
rz(-pi) q[2];
rz(2.4268552) q[3];
sx q[3];
rz(-0.5849896) q[3];
sx q[3];
rz(1.8470256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-2.5626101) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-1.0639169) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.3508505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3063072) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(-1.6552734) q[0];
x q[1];
rz(0.66019085) q[2];
sx q[2];
rz(-1.635951) q[2];
sx q[2];
rz(0.24812631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1119712) q[1];
sx q[1];
rz(-2.0051149) q[1];
sx q[1];
rz(-3.1236468) q[1];
rz(-pi) q[2];
rz(2.001858) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.110048) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(0.44719493) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34925941) q[0];
sx q[0];
rz(-2.0418641) q[0];
sx q[0];
rz(0.40339289) q[0];
x q[1];
rz(-2.7197003) q[2];
sx q[2];
rz(-1.5932398) q[2];
sx q[2];
rz(-2.9107712) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.604653) q[1];
sx q[1];
rz(-2.5250658) q[1];
sx q[1];
rz(-1.489868) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-0.51207536) q[1];
sx q[1];
rz(-2.8072186) q[1];
sx q[1];
rz(1.0288815) q[1];
rz(1.8722664) q[2];
sx q[2];
rz(-1.3133247) q[2];
sx q[2];
rz(1.7016344) q[2];
rz(-0.11733304) q[3];
sx q[3];
rz(-1.7775848) q[3];
sx q[3];
rz(2.1669273) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
