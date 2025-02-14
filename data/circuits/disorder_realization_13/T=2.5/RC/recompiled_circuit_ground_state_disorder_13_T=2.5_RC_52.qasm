OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(1.2989651) q[0];
rz(2.9453912) q[1];
sx q[1];
rz(-1.8008404) q[1];
sx q[1];
rz(0.10032108) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8232788) q[0];
sx q[0];
rz(-2.4779076) q[0];
sx q[0];
rz(0.084245988) q[0];
rz(0.47122987) q[2];
sx q[2];
rz(-1.2696075) q[2];
sx q[2];
rz(2.9848841) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92857526) q[1];
sx q[1];
rz(-1.2963821) q[1];
sx q[1];
rz(-1.3540512) q[1];
rz(-pi) q[2];
rz(-2.6907608) q[3];
sx q[3];
rz(-0.63967645) q[3];
sx q[3];
rz(-2.0995288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4002865) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(2.3347704) q[2];
rz(0.72871366) q[3];
sx q[3];
rz(-2.3829134) q[3];
sx q[3];
rz(-1.9369283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5173986) q[0];
sx q[0];
rz(-2.875858) q[0];
sx q[0];
rz(-0.30701315) q[0];
rz(1.864805) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(-1.101864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1043423) q[0];
sx q[0];
rz(-0.58339632) q[0];
sx q[0];
rz(1.8895545) q[0];
x q[1];
rz(0.33874933) q[2];
sx q[2];
rz(-2.951606) q[2];
sx q[2];
rz(2.7471586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3668573) q[1];
sx q[1];
rz(-0.22379074) q[1];
sx q[1];
rz(-0.91694497) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83294932) q[3];
sx q[3];
rz(-2.2112339) q[3];
sx q[3];
rz(-2.4082977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.028194204) q[2];
sx q[2];
rz(-2.5598309) q[2];
sx q[2];
rz(-3.0451575) q[2];
rz(-2.9891369) q[3];
sx q[3];
rz(-1.6343445) q[3];
sx q[3];
rz(2.4436387) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0518799) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-0.40019792) q[0];
rz(-1.5125037) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(0.2581183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5921051) q[0];
sx q[0];
rz(-1.4319812) q[0];
sx q[0];
rz(-1.8068683) q[0];
x q[1];
rz(0.099868593) q[2];
sx q[2];
rz(-1.6503008) q[2];
sx q[2];
rz(0.85209633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3052747) q[1];
sx q[1];
rz(-1.5757676) q[1];
sx q[1];
rz(0.012955091) q[1];
rz(-0.6155562) q[3];
sx q[3];
rz(-0.29272348) q[3];
sx q[3];
rz(-2.9396201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.021598024) q[2];
sx q[2];
rz(-1.1979411) q[2];
sx q[2];
rz(3.1094587) q[2];
rz(2.8739127) q[3];
sx q[3];
rz(-1.6240424) q[3];
sx q[3];
rz(1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5716008) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-3.0573523) q[0];
rz(0.035942297) q[1];
sx q[1];
rz(-3.1075931) q[1];
sx q[1];
rz(-0.34119225) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81927903) q[0];
sx q[0];
rz(-1.8321165) q[0];
sx q[0];
rz(-0.65066353) q[0];
rz(-pi) q[1];
rz(1.313339) q[2];
sx q[2];
rz(-0.99960589) q[2];
sx q[2];
rz(-1.5866304) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2349639) q[1];
sx q[1];
rz(-0.90032691) q[1];
sx q[1];
rz(-0.57256289) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61670035) q[3];
sx q[3];
rz(-2.4444207) q[3];
sx q[3];
rz(2.8215849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42698947) q[2];
sx q[2];
rz(-2.0689071) q[2];
sx q[2];
rz(1.6061456) q[2];
rz(2.8793907) q[3];
sx q[3];
rz(-1.5235135) q[3];
sx q[3];
rz(0.95482993) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75761211) q[0];
sx q[0];
rz(-0.42344991) q[0];
sx q[0];
rz(2.0641932) q[0];
rz(2.7491838) q[1];
sx q[1];
rz(-3.0632186) q[1];
sx q[1];
rz(2.1108625) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.01973132) q[0];
sx q[0];
rz(-1.1543589) q[0];
sx q[0];
rz(0.58953489) q[0];
rz(-pi) q[1];
rz(-0.18056492) q[2];
sx q[2];
rz(-2.0624277) q[2];
sx q[2];
rz(-2.5485843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2490152) q[1];
sx q[1];
rz(-1.4059781) q[1];
sx q[1];
rz(1.3457764) q[1];
x q[2];
rz(-2.137759) q[3];
sx q[3];
rz(-0.58708411) q[3];
sx q[3];
rz(2.7239012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3173759) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(0.83924323) q[2];
rz(-0.9134891) q[3];
sx q[3];
rz(-1.313611) q[3];
sx q[3];
rz(-2.998108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4588673) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(-1.4759901) q[0];
rz(-2.7492375) q[1];
sx q[1];
rz(-2.0456435) q[1];
sx q[1];
rz(2.5501693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2570193) q[0];
sx q[0];
rz(-1.2019751) q[0];
sx q[0];
rz(-0.9414282) q[0];
rz(-pi) q[1];
rz(2.8663581) q[2];
sx q[2];
rz(-1.7448336) q[2];
sx q[2];
rz(-1.7013719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6900031) q[1];
sx q[1];
rz(-1.065548) q[1];
sx q[1];
rz(-3.0662159) q[1];
x q[2];
rz(-0.87402451) q[3];
sx q[3];
rz(-1.2748147) q[3];
sx q[3];
rz(-3.0324453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8412987) q[2];
sx q[2];
rz(-2.4755307) q[2];
sx q[2];
rz(-0.57233125) q[2];
rz(2.9193997) q[3];
sx q[3];
rz(-0.43155813) q[3];
sx q[3];
rz(0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.38729024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(2.733316) q[0];
rz(0.73221842) q[1];
sx q[1];
rz(-0.1258985) q[1];
sx q[1];
rz(2.8439723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10935303) q[0];
sx q[0];
rz(-1.223757) q[0];
sx q[0];
rz(0.45022398) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6082798) q[2];
sx q[2];
rz(-2.6465694) q[2];
sx q[2];
rz(-1.9062454) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0113161) q[1];
sx q[1];
rz(-1.374048) q[1];
sx q[1];
rz(1.9185683) q[1];
rz(-pi) q[2];
rz(2.9651661) q[3];
sx q[3];
rz(-0.71143141) q[3];
sx q[3];
rz(-2.3030858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.97062651) q[2];
sx q[2];
rz(-1.2622702) q[2];
sx q[2];
rz(-0.69620281) q[2];
rz(2.2512186) q[3];
sx q[3];
rz(-1.1768769) q[3];
sx q[3];
rz(1.4674998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.17906976) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(-0.21275511) q[0];
rz(2.6720324) q[1];
sx q[1];
rz(-0.9649562) q[1];
sx q[1];
rz(2.3874217) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6109723) q[0];
sx q[0];
rz(-1.3830796) q[0];
sx q[0];
rz(-1.8646445) q[0];
rz(0.55108549) q[2];
sx q[2];
rz(-2.7646825) q[2];
sx q[2];
rz(-2.2152679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8629145) q[1];
sx q[1];
rz(-0.92611852) q[1];
sx q[1];
rz(1.7134604) q[1];
rz(0.92408224) q[3];
sx q[3];
rz(-1.0093186) q[3];
sx q[3];
rz(1.1738861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.135005) q[2];
sx q[2];
rz(-0.90234119) q[2];
sx q[2];
rz(-0.78224409) q[2];
rz(-1.4462645) q[3];
sx q[3];
rz(-2.5952314) q[3];
sx q[3];
rz(-2.7915891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0924454) q[0];
sx q[0];
rz(-2.6706084) q[0];
sx q[0];
rz(-2.1771722) q[0];
rz(-1.8655221) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.6395578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1196331) q[0];
sx q[0];
rz(-0.30888452) q[0];
sx q[0];
rz(-0.85794411) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1310857) q[2];
sx q[2];
rz(-0.88958101) q[2];
sx q[2];
rz(-1.1352254) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.137845) q[1];
sx q[1];
rz(-1.5478351) q[1];
sx q[1];
rz(-2.0883043) q[1];
rz(0.067599452) q[3];
sx q[3];
rz(-1.7177594) q[3];
sx q[3];
rz(-2.802668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2532578) q[2];
sx q[2];
rz(-1.2829245) q[2];
sx q[2];
rz(0.9453195) q[2];
rz(2.3156598) q[3];
sx q[3];
rz(-1.632894) q[3];
sx q[3];
rz(0.53502423) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0655521) q[0];
sx q[0];
rz(-1.1689508) q[0];
sx q[0];
rz(-2.3384576) q[0];
rz(-1.5638634) q[1];
sx q[1];
rz(-1.4811265) q[1];
sx q[1];
rz(-0.28958431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3428925) q[0];
sx q[0];
rz(-3.0353167) q[0];
sx q[0];
rz(-1.6864087) q[0];
rz(-pi) q[1];
rz(-1.8990191) q[2];
sx q[2];
rz(-1.8716836) q[2];
sx q[2];
rz(-0.93133486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9152569) q[1];
sx q[1];
rz(-0.42789547) q[1];
sx q[1];
rz(-1.0029327) q[1];
x q[2];
rz(-2.2749994) q[3];
sx q[3];
rz(-2.7150053) q[3];
sx q[3];
rz(-2.085919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.76558602) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(0.96013367) q[2];
rz(-0.58297408) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(-0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.690602) q[0];
sx q[0];
rz(-1.4324181) q[0];
sx q[0];
rz(1.6313534) q[0];
rz(-0.040738978) q[1];
sx q[1];
rz(-2.4650885) q[1];
sx q[1];
rz(-3.0104641) q[1];
rz(-2.8340451) q[2];
sx q[2];
rz(-1.1196305) q[2];
sx q[2];
rz(1.011871) q[2];
rz(1.4996281) q[3];
sx q[3];
rz(-1.6215848) q[3];
sx q[3];
rz(3.0231089) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
