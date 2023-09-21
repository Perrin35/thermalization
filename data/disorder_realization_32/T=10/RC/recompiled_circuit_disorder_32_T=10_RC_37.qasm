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
rz(-1.0892692) q[0];
sx q[0];
rz(-2.9805592) q[0];
rz(1.610202) q[1];
sx q[1];
rz(2.664497) q[1];
sx q[1];
rz(8.9283979) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3408605) q[0];
sx q[0];
rz(-0.98156089) q[0];
sx q[0];
rz(0.13829921) q[0];
x q[1];
rz(-1.7726937) q[2];
sx q[2];
rz(-1.9041833) q[2];
sx q[2];
rz(2.8436529) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15064083) q[1];
sx q[1];
rz(-2.4283263) q[1];
sx q[1];
rz(0.13891061) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6707889) q[3];
sx q[3];
rz(-1.1649719) q[3];
sx q[3];
rz(1.0909181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6266142) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(-1.3927762) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(-1.1028517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9614812) q[0];
sx q[0];
rz(-2.8948088) q[0];
sx q[0];
rz(-2.6632975) q[0];
rz(-pi) q[1];
rz(-2.8385542) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(2.9532202) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9175375) q[1];
sx q[1];
rz(-1.1174669) q[1];
sx q[1];
rz(-0.097150306) q[1];
rz(-0.79264499) q[3];
sx q[3];
rz(-2.9229197) q[3];
sx q[3];
rz(-1.6817026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(-2.0593026) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796233) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(-2.9586155) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.9972237) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367033) q[0];
sx q[0];
rz(-1.0495782) q[0];
sx q[0];
rz(0.46023603) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7230354) q[2];
sx q[2];
rz(-1.7757799) q[2];
sx q[2];
rz(0.95450729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.548) q[1];
sx q[1];
rz(-1.9152194) q[1];
sx q[1];
rz(0.48732948) q[1];
x q[2];
rz(-2.8539682) q[3];
sx q[3];
rz(-1.5737185) q[3];
sx q[3];
rz(0.50141108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(0.90467492) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-2.3864746) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6782137) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-2.082086) q[1];
sx q[1];
rz(-0.99266565) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794605) q[0];
sx q[0];
rz(-1.3068763) q[0];
sx q[0];
rz(-2.604565) q[0];
rz(-pi) q[1];
rz(0.56447345) q[2];
sx q[2];
rz(-2.1285004) q[2];
sx q[2];
rz(-1.0587453) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.48651931) q[1];
sx q[1];
rz(-1.4078948) q[1];
sx q[1];
rz(-1.4763114) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1487001) q[3];
sx q[3];
rz(-1.9979457) q[3];
sx q[3];
rz(-0.072673365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59584004) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(-1.9013566) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-0.3516745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2857916) q[0];
sx q[0];
rz(-1.283678) q[0];
sx q[0];
rz(2.48824) q[0];
rz(-pi) q[1];
rz(-0.50260966) q[2];
sx q[2];
rz(-2.1447499) q[2];
sx q[2];
rz(0.9290907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5366718) q[1];
sx q[1];
rz(-2.7465638) q[1];
sx q[1];
rz(3.1035963) q[1];
rz(-pi) q[2];
rz(-0.3211) q[3];
sx q[3];
rz(-1.0874815) q[3];
sx q[3];
rz(1.4534284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(2.5999516) q[0];
rz(0.60846865) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(-2.0419962) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0352286) q[0];
sx q[0];
rz(-1.117825) q[0];
sx q[0];
rz(-0.41082541) q[0];
x q[1];
rz(-1.8855368) q[2];
sx q[2];
rz(-0.55329269) q[2];
sx q[2];
rz(1.7994583) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82367831) q[1];
sx q[1];
rz(-0.70033973) q[1];
sx q[1];
rz(-0.27842303) q[1];
rz(-2.3533456) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(-1.6227674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5114484) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0866078) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(-0.64487547) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(0.15730102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511714) q[0];
sx q[0];
rz(-0.84050814) q[0];
sx q[0];
rz(1.4448326) q[0];
rz(-pi) q[1];
rz(-2.925161) q[2];
sx q[2];
rz(-0.84605233) q[2];
sx q[2];
rz(-0.12491465) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0485222) q[1];
sx q[1];
rz(-1.7659566) q[1];
sx q[1];
rz(2.1919769) q[1];
rz(-pi) q[2];
rz(1.2792475) q[3];
sx q[3];
rz(-0.96372094) q[3];
sx q[3];
rz(0.32015043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(-2.9390745) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4956932) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(2.8483134) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5787443) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(-1.5777274) q[0];
x q[1];
rz(2.9759334) q[2];
sx q[2];
rz(-2.207901) q[2];
sx q[2];
rz(-3.0050302) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5756702) q[1];
sx q[1];
rz(-0.88469632) q[1];
sx q[1];
rz(-2.5887606) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(-2.8914715) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(-2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0042689) q[0];
sx q[0];
rz(-0.22432888) q[0];
sx q[0];
rz(-1.0639169) q[0];
rz(-2.6990199) q[1];
sx q[1];
rz(-1.7174218) q[1];
sx q[1];
rz(-1.7907422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8838053) q[0];
sx q[0];
rz(-1.4865849) q[0];
sx q[0];
rz(3.0620831) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6532134) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(-1.8694307) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0693448) q[1];
sx q[1];
rz(-2.7069271) q[1];
sx q[1];
rz(1.5321295) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1397347) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(-2.7800625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-2.6943977) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31496012) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(-2.3610624) q[0];
rz(-2.7138117) q[1];
sx q[1];
rz(-2.784274) q[1];
sx q[1];
rz(-2.9719877) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0376301) q[0];
sx q[0];
rz(-2.5314405) q[0];
sx q[0];
rz(2.2274341) q[0];
x q[1];
rz(0.42189235) q[2];
sx q[2];
rz(-1.5932398) q[2];
sx q[2];
rz(0.23082146) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1093724) q[1];
sx q[1];
rz(-1.5240372) q[1];
sx q[1];
rz(-2.1857775) q[1];
x q[2];
rz(0.1428991) q[3];
sx q[3];
rz(-1.2522962) q[3];
sx q[3];
rz(-2.3058476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85049373) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(-2.4492241) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904084) q[0];
sx q[0];
rz(-1.6155227) q[0];
sx q[0];
rz(0.68646705) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(1.2693263) q[2];
sx q[2];
rz(-1.8282679) q[2];
sx q[2];
rz(-1.4399583) q[2];
rz(1.3626171) q[3];
sx q[3];
rz(-1.4559742) q[3];
sx q[3];
rz(-2.5696587) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];