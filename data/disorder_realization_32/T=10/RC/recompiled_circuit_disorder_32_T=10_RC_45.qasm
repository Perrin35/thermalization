OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.035064) q[0];
sx q[0];
rz(-2.0523235) q[0];
sx q[0];
rz(2.9805592) q[0];
rz(-1.5313907) q[1];
sx q[1];
rz(-2.664497) q[1];
sx q[1];
rz(2.6452126) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0954477) q[0];
sx q[0];
rz(-0.60337043) q[0];
sx q[0];
rz(-1.3674111) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6167294) q[2];
sx q[2];
rz(-0.3877936) q[2];
sx q[2];
rz(-0.85688574) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6160994) q[1];
sx q[1];
rz(-1.4800737) q[1];
sx q[1];
rz(2.4331122) q[1];
rz(1.6707889) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(-2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51497841) q[2];
sx q[2];
rz(-2.2815621) q[2];
sx q[2];
rz(2.853945) q[2];
rz(-1.3927762) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(2.0387409) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2709133) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(2.5640008) q[0];
rz(-1.5354935) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(-1.5637406) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751511) q[0];
sx q[0];
rz(-1.4581212) q[0];
sx q[0];
rz(-2.9215647) q[0];
x q[1];
rz(2.4096476) q[2];
sx q[2];
rz(-0.39790301) q[2];
sx q[2];
rz(1.0674455) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9175375) q[1];
sx q[1];
rz(-2.0241258) q[1];
sx q[1];
rz(-0.097150306) q[1];
rz(-pi) q[2];
rz(1.413826) q[3];
sx q[3];
rz(-1.7236712) q[3];
sx q[3];
rz(2.4863941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-2.1652174) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(-1.9783431) q[3];
sx q[3];
rz(-1.9415559) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.0619693) q[0];
sx q[0];
rz(-0.29456961) q[0];
sx q[0];
rz(-2.9586155) q[0];
rz(0.043116365) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(-1.9972237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9632918) q[0];
sx q[0];
rz(-2.460647) q[0];
sx q[0];
rz(0.91239022) q[0];
rz(-0.41855721) q[2];
sx q[2];
rz(-1.7757799) q[2];
sx q[2];
rz(0.95450729) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.548) q[1];
sx q[1];
rz(-1.2263733) q[1];
sx q[1];
rz(-0.48732948) q[1];
x q[2];
rz(0.010300962) q[3];
sx q[3];
rz(-0.28763887) q[3];
sx q[3];
rz(-2.0623296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(2.2369177) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-0.42680877) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-0.97020522) q[0];
sx q[0];
rz(0.71587193) q[0];
rz(0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(0.99266565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0794605) q[0];
sx q[0];
rz(-1.8347164) q[0];
sx q[0];
rz(2.604565) q[0];
rz(2.2797212) q[2];
sx q[2];
rz(-0.77152354) q[2];
sx q[2];
rz(2.9574403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0992972) q[1];
sx q[1];
rz(-2.9534833) q[1];
sx q[1];
rz(-0.52109615) q[1];
rz(-2.643814) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(-1.9073515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0830393) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.7129718) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.5139791) q[3];
sx q[3];
rz(1.9184453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(2.3748421) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(0.3516745) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.283678) q[0];
sx q[0];
rz(0.65335269) q[0];
rz(-0.50260966) q[2];
sx q[2];
rz(-2.1447499) q[2];
sx q[2];
rz(-2.212502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5366718) q[1];
sx q[1];
rz(-0.39502883) q[1];
sx q[1];
rz(-0.037996304) q[1];
x q[2];
rz(1.0655754) q[3];
sx q[3];
rz(-1.8540283) q[3];
sx q[3];
rz(0.035988228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19796431) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(-0.42896459) q[3];
sx q[3];
rz(-0.52754378) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-2.5999516) q[0];
rz(2.533124) q[1];
sx q[1];
rz(-0.21921961) q[1];
sx q[1];
rz(2.0419962) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8655411) q[0];
sx q[0];
rz(-1.2035032) q[0];
sx q[0];
rz(-1.0827351) q[0];
rz(-1.0397644) q[2];
sx q[2];
rz(-1.4073939) q[2];
sx q[2];
rz(2.642717) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9601599) q[1];
sx q[1];
rz(-2.2391041) q[1];
sx q[1];
rz(-1.3431576) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8058365) q[3];
sx q[3];
rz(-0.81695518) q[3];
sx q[3];
rz(0.18248617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.6116066) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-2.6480643) q[1];
sx q[1];
rz(0.15730102) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8511714) q[0];
sx q[0];
rz(-2.3010845) q[0];
sx q[0];
rz(1.4448326) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3072725) q[2];
sx q[2];
rz(-1.7322707) q[2];
sx q[2];
rz(1.8404567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.09307043) q[1];
sx q[1];
rz(-1.375636) q[1];
sx q[1];
rz(-2.1919769) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5141684) q[3];
sx q[3];
rz(-1.809123) q[3];
sx q[3];
rz(2.0605007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9592231) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(3.0977141) q[3];
sx q[3];
rz(-1.9079804) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(-3.1307401) q[0];
rz(0.27663484) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(0.29327926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56284833) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(1.5638652) q[0];
rz(1.3515527) q[2];
sx q[2];
rz(-2.4862137) q[2];
sx q[2];
rz(-3.0041681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7640904) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(2.3368895) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46383143) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(-2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(-0.088767178) q[2];
rz(2.8914715) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(-2.5626101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
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
rz(1.7907422) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1254743) q[0];
sx q[0];
rz(-3.0258412) q[0];
sx q[0];
rz(-0.8158169) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0355989) q[2];
sx q[2];
rz(-0.66291891) q[2];
sx q[2];
rz(1.4063327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1119712) q[1];
sx q[1];
rz(-1.1364778) q[1];
sx q[1];
rz(0.017945826) q[1];
rz(1.1397347) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(-2.7800625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.110048) q[2];
sx q[2];
rz(-1.2908547) q[2];
sx q[2];
rz(-0.44719493) q[2];
rz(-1.3859008) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(-2.6269004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(-0.78053027) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0302089) q[0];
sx q[0];
rz(-1.2134524) q[0];
sx q[0];
rz(2.0765199) q[0];
rz(-pi) q[1];
rz(-0.054762997) q[2];
sx q[2];
rz(-2.7191396) q[2];
sx q[2];
rz(-1.7516608) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5055713) q[1];
sx q[1];
rz(-0.95658703) q[1];
sx q[1];
rz(0.057227055) q[1];
rz(-pi) q[2];
rz(-1.1630837) q[3];
sx q[3];
rz(-2.793503) q[3];
sx q[3];
rz(0.40504328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.1824181) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(-0.46323562) q[3];
sx q[3];
rz(-1.4866456) q[3];
sx q[3];
rz(-2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2904084) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(-2.6295173) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(2.8725273) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
rz(-2.0797707) q[3];
sx q[3];
rz(-0.23734262) q[3];
sx q[3];
rz(1.6457641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
