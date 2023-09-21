OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(-0.0097302516) q[1];
sx q[1];
rz(-1.4571804) q[1];
sx q[1];
rz(1.943346) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449013) q[0];
sx q[0];
rz(-2.9400819) q[0];
sx q[0];
rz(-2.6411396) q[0];
rz(2.8134402) q[2];
sx q[2];
rz(-2.3386764) q[2];
sx q[2];
rz(-0.83628718) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4858422) q[1];
sx q[1];
rz(-1.2409004) q[1];
sx q[1];
rz(2.5962023) q[1];
rz(0.046273307) q[3];
sx q[3];
rz(-1.4001906) q[3];
sx q[3];
rz(1.6149278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1951695) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(0.18134376) q[2];
rz(2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-2.3852824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392035) q[0];
sx q[0];
rz(-2.8518682) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(-0.50239262) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(-1.5418672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4138448) q[0];
sx q[0];
rz(-1.9808597) q[0];
sx q[0];
rz(-0.64322612) q[0];
x q[1];
rz(0.16427152) q[2];
sx q[2];
rz(-1.1195682) q[2];
sx q[2];
rz(-0.8115561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.220037) q[1];
sx q[1];
rz(-1.5350635) q[1];
sx q[1];
rz(1.5807371) q[1];
rz(0.48331355) q[3];
sx q[3];
rz(-1.1477074) q[3];
sx q[3];
rz(1.9139293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(-0.85033068) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.458228) q[3];
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
rz(2.8787815) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(0.61022726) q[0];
rz(-1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-0.99197018) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6863166) q[0];
sx q[0];
rz(-2.1222097) q[0];
sx q[0];
rz(-0.21689143) q[0];
rz(-pi) q[1];
x q[1];
rz(2.473258) q[2];
sx q[2];
rz(-1.4112345) q[2];
sx q[2];
rz(2.4192686) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5132644) q[1];
sx q[1];
rz(-0.4497512) q[1];
sx q[1];
rz(0.20723923) q[1];
rz(-2.6716258) q[3];
sx q[3];
rz(-2.6372058) q[3];
sx q[3];
rz(-0.72987635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(2.8733011) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-0.18850732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2878993) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(-2.4515117) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(-2.4437723) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9245802) q[0];
sx q[0];
rz(-0.092886535) q[0];
sx q[0];
rz(-1.8266982) q[0];
rz(-pi) q[1];
rz(-1.1890829) q[2];
sx q[2];
rz(-2.082798) q[2];
sx q[2];
rz(-2.2005759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88155014) q[1];
sx q[1];
rz(-1.3870194) q[1];
sx q[1];
rz(1.0764513) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29699765) q[3];
sx q[3];
rz(-1.9747435) q[3];
sx q[3];
rz(-1.7904074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(-1.6323803) q[2];
rz(-0.40431067) q[3];
sx q[3];
rz(-2.4590838) q[3];
sx q[3];
rz(1.4782762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8835835) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(-1.0255381) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(-0.62932032) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.660491) q[0];
sx q[0];
rz(-1.2718624) q[0];
sx q[0];
rz(2.7116508) q[0];
rz(-0.25482486) q[2];
sx q[2];
rz(-0.91081589) q[2];
sx q[2];
rz(-0.13435907) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3931261) q[1];
sx q[1];
rz(-1.6854291) q[1];
sx q[1];
rz(2.6266891) q[1];
rz(-pi) q[2];
rz(-2.1719645) q[3];
sx q[3];
rz(-1.789635) q[3];
sx q[3];
rz(0.18274433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59262529) q[2];
sx q[2];
rz(-0.36965814) q[2];
sx q[2];
rz(0.34234753) q[2];
rz(1.3458378) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(-0.19601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.65258566) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(2.956399) q[0];
rz(1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(1.7746183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.470984) q[0];
sx q[0];
rz(-2.1489722) q[0];
sx q[0];
rz(-1.3060119) q[0];
rz(-pi) q[1];
rz(-0.30250678) q[2];
sx q[2];
rz(-2.2777646) q[2];
sx q[2];
rz(-0.44755852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1176599) q[1];
sx q[1];
rz(-1.6624833) q[1];
sx q[1];
rz(1.7232643) q[1];
rz(-pi) q[2];
rz(-2.572445) q[3];
sx q[3];
rz(-1.2043118) q[3];
sx q[3];
rz(-1.5018499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8967445) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(0.81378716) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8900523) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.3758804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863231) q[0];
sx q[0];
rz(-0.77676847) q[0];
sx q[0];
rz(-1.1458678) q[0];
rz(1.6872348) q[2];
sx q[2];
rz(-2.6652626) q[2];
sx q[2];
rz(-1.6778698) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3657276) q[1];
sx q[1];
rz(-2.9850246) q[1];
sx q[1];
rz(-0.74128976) q[1];
rz(-pi) q[2];
rz(1.901058) q[3];
sx q[3];
rz(-1.6263279) q[3];
sx q[3];
rz(0.33952573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4574796) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(2.4105371) q[2];
rz(0.11080065) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(-0.69537648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59657997) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(-0.73356432) q[0];
rz(-0.60797524) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9600784) q[0];
sx q[0];
rz(-2.5900822) q[0];
sx q[0];
rz(-2.0803204) q[0];
rz(-pi) q[1];
rz(2.4085326) q[2];
sx q[2];
rz(-0.86793938) q[2];
sx q[2];
rz(-1.5066063) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9573606) q[1];
sx q[1];
rz(-2.3174274) q[1];
sx q[1];
rz(0.25914534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0844564) q[3];
sx q[3];
rz(-2.497695) q[3];
sx q[3];
rz(0.38734303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.3092821) q[2];
sx q[2];
rz(-1.7162494) q[2];
rz(1.6783293) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5995246) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(-1.19338) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.3648938) q[1];
sx q[1];
rz(-2.1967922) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9076618) q[0];
sx q[0];
rz(-0.91066563) q[0];
sx q[0];
rz(-1.075676) q[0];
x q[1];
rz(0.94061942) q[2];
sx q[2];
rz(-0.47626469) q[2];
sx q[2];
rz(0.95048743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8972842) q[1];
sx q[1];
rz(-2.8123887) q[1];
sx q[1];
rz(2.5876178) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9267843) q[3];
sx q[3];
rz(-1.3857406) q[3];
sx q[3];
rz(-1.9389648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.21645674) q[2];
sx q[2];
rz(-1.0711121) q[2];
sx q[2];
rz(0.28820583) q[2];
rz(0.47973412) q[3];
sx q[3];
rz(-1.0498485) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11186803) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(-2.2286041) q[0];
rz(0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(-2.250681) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4961632) q[0];
sx q[0];
rz(-1.9053639) q[0];
sx q[0];
rz(-1.8339001) q[0];
rz(-pi) q[1];
rz(2.8577523) q[2];
sx q[2];
rz(-0.97728697) q[2];
sx q[2];
rz(2.3317091) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18894698) q[1];
sx q[1];
rz(-1.2332321) q[1];
sx q[1];
rz(-1.7811716) q[1];
rz(-pi) q[2];
rz(2.7632305) q[3];
sx q[3];
rz(-0.20212691) q[3];
sx q[3];
rz(-2.9581021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23641071) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(0.50160828) q[2];
rz(1.2891399) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63507737) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(1.1322017) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-2.6272527) q[2];
sx q[2];
rz(-2.8375576) q[2];
sx q[2];
rz(0.67630771) q[2];
rz(0.049384762) q[3];
sx q[3];
rz(-1.0247083) q[3];
sx q[3];
rz(0.68030737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];