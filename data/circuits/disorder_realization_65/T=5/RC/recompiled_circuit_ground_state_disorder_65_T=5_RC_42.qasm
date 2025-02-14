OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.68706694) q[0];
sx q[0];
rz(-1.878976) q[0];
sx q[0];
rz(2.4071121) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(-1.2999111) q[1];
sx q[1];
rz(-2.1631961) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1837195) q[0];
sx q[0];
rz(-1.5629725) q[0];
sx q[0];
rz(1.70723) q[0];
rz(-2.4327203) q[2];
sx q[2];
rz(-1.1507431) q[2];
sx q[2];
rz(1.5576897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96127993) q[1];
sx q[1];
rz(-2.790457) q[1];
sx q[1];
rz(1.9957955) q[1];
rz(-pi) q[2];
rz(2.0693073) q[3];
sx q[3];
rz(-1.5007334) q[3];
sx q[3];
rz(0.69341502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3413099) q[2];
sx q[2];
rz(-1.492123) q[2];
sx q[2];
rz(0.71551234) q[2];
rz(-1.4095151) q[3];
sx q[3];
rz(-2.2655497) q[3];
sx q[3];
rz(1.5040262) q[3];
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
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3187123) q[0];
sx q[0];
rz(-2.5647698) q[0];
sx q[0];
rz(-3.1067644) q[0];
rz(0.71827978) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(1.1911596) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1930192) q[0];
sx q[0];
rz(-0.64136845) q[0];
sx q[0];
rz(-1.3032342) q[0];
rz(2.4592752) q[2];
sx q[2];
rz(-1.6860262) q[2];
sx q[2];
rz(-0.29733411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71514788) q[1];
sx q[1];
rz(-0.99438018) q[1];
sx q[1];
rz(-1.9020686) q[1];
x q[2];
rz(-2.0091811) q[3];
sx q[3];
rz(-1.7436641) q[3];
sx q[3];
rz(0.57314429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3388153) q[2];
sx q[2];
rz(-1.2884527) q[2];
sx q[2];
rz(-0.046099376) q[2];
rz(0.80373803) q[3];
sx q[3];
rz(-1.2396783) q[3];
sx q[3];
rz(0.55155915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4285202) q[0];
sx q[0];
rz(-1.5794733) q[0];
sx q[0];
rz(-0.61306104) q[0];
rz(-2.8763981) q[1];
sx q[1];
rz(-1.801633) q[1];
sx q[1];
rz(3.0934966) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4713132) q[0];
sx q[0];
rz(-1.0889772) q[0];
sx q[0];
rz(3.1387718) q[0];
rz(0.71718054) q[2];
sx q[2];
rz(-2.4240026) q[2];
sx q[2];
rz(-2.8924475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.887943) q[1];
sx q[1];
rz(-3.0620975) q[1];
sx q[1];
rz(-0.84913932) q[1];
rz(-1.223391) q[3];
sx q[3];
rz(-1.9022397) q[3];
sx q[3];
rz(-3.0104266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4464104) q[2];
sx q[2];
rz(-1.9789275) q[2];
sx q[2];
rz(2.7282558) q[2];
rz(-0.6684331) q[3];
sx q[3];
rz(-1.8139402) q[3];
sx q[3];
rz(1.7641164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.9628825) q[0];
sx q[0];
rz(-0.2186192) q[0];
sx q[0];
rz(2.9982153) q[0];
rz(-2.7168221) q[1];
sx q[1];
rz(-1.9353341) q[1];
sx q[1];
rz(0.43407789) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33725564) q[0];
sx q[0];
rz(-1.6640264) q[0];
sx q[0];
rz(-3.1033959) q[0];
x q[1];
rz(0.6521449) q[2];
sx q[2];
rz(-1.6639478) q[2];
sx q[2];
rz(0.31983122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0632759) q[1];
sx q[1];
rz(-2.6022291) q[1];
sx q[1];
rz(-1.2961948) q[1];
x q[2];
rz(-2.3301773) q[3];
sx q[3];
rz(-1.7219117) q[3];
sx q[3];
rz(-2.9600157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8726337) q[2];
sx q[2];
rz(-0.7577529) q[2];
sx q[2];
rz(2.7244275) q[2];
rz(0.68552351) q[3];
sx q[3];
rz(-1.439582) q[3];
sx q[3];
rz(-3.0912257) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79614574) q[0];
sx q[0];
rz(-2.8247927) q[0];
sx q[0];
rz(-1.6337122) q[0];
rz(2.4729074) q[1];
sx q[1];
rz(-1.1983127) q[1];
sx q[1];
rz(2.0100458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.749866) q[0];
sx q[0];
rz(-1.6035115) q[0];
sx q[0];
rz(-2.1689438) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0102398) q[2];
sx q[2];
rz(-1.2850637) q[2];
sx q[2];
rz(3.07522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9100489) q[1];
sx q[1];
rz(-1.1041512) q[1];
sx q[1];
rz(-1.3620767) q[1];
rz(-0.82238166) q[3];
sx q[3];
rz(-0.17381343) q[3];
sx q[3];
rz(-1.8552775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23182997) q[2];
sx q[2];
rz(-2.6015687) q[2];
sx q[2];
rz(-2.6798999) q[2];
rz(-0.26990226) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(-2.4902978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57399026) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(-1.061576) q[0];
rz(-2.4827982) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(-0.40491358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39570582) q[0];
sx q[0];
rz(-2.2721724) q[0];
sx q[0];
rz(1.7084136) q[0];
rz(-1.986386) q[2];
sx q[2];
rz(-1.0672369) q[2];
sx q[2];
rz(-1.0114618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3337219) q[1];
sx q[1];
rz(-2.0757339) q[1];
sx q[1];
rz(-2.9823279) q[1];
rz(-1.177321) q[3];
sx q[3];
rz(-0.67613908) q[3];
sx q[3];
rz(-2.2666988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0460661) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(-1.3859762) q[2];
rz(-2.4447377) q[3];
sx q[3];
rz(-0.98074061) q[3];
sx q[3];
rz(0.81010404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(2.5974469) q[0];
sx q[0];
rz(-2.5283041) q[0];
sx q[0];
rz(-2.5469653) q[0];
rz(-2.0095297) q[1];
sx q[1];
rz(-1.7375528) q[1];
sx q[1];
rz(2.6304257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5385732) q[0];
sx q[0];
rz(-1.5165189) q[0];
sx q[0];
rz(0.57247573) q[0];
rz(-pi) q[1];
rz(2.2726502) q[2];
sx q[2];
rz(-2.7232183) q[2];
sx q[2];
rz(-1.4792031) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6936924) q[1];
sx q[1];
rz(-1.3096454) q[1];
sx q[1];
rz(3.0413701) q[1];
x q[2];
rz(2.0603544) q[3];
sx q[3];
rz(-0.900002) q[3];
sx q[3];
rz(-1.9157992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96659294) q[2];
sx q[2];
rz(-2.9006557) q[2];
sx q[2];
rz(2.8540376) q[2];
rz(2.0368841) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(0.12565676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.48968807) q[0];
sx q[0];
rz(-0.37864417) q[0];
sx q[0];
rz(2.4878159) q[0];
rz(-2.6405624) q[1];
sx q[1];
rz(-0.72622314) q[1];
sx q[1];
rz(0.78530606) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619072) q[0];
sx q[0];
rz(-2.5568049) q[0];
sx q[0];
rz(-0.25081046) q[0];
rz(1.1210322) q[2];
sx q[2];
rz(-1.0503979) q[2];
sx q[2];
rz(0.880151) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.056902) q[1];
sx q[1];
rz(-2.4246019) q[1];
sx q[1];
rz(-0.3221953) q[1];
x q[2];
rz(1.4887962) q[3];
sx q[3];
rz(-1.3424216) q[3];
sx q[3];
rz(-2.5395951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.952897) q[2];
sx q[2];
rz(-1.5263824) q[2];
sx q[2];
rz(2.7719899) q[2];
rz(2.8772307) q[3];
sx q[3];
rz(-1.8601469) q[3];
sx q[3];
rz(3.1407691) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13067506) q[0];
sx q[0];
rz(-0.70931804) q[0];
sx q[0];
rz(1.6226907) q[0];
rz(1.2394637) q[1];
sx q[1];
rz(-1.112097) q[1];
sx q[1];
rz(-0.94737238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0305875) q[0];
sx q[0];
rz(-1.9698785) q[0];
sx q[0];
rz(-3.0754979) q[0];
rz(-pi) q[1];
rz(0.36115335) q[2];
sx q[2];
rz(-1.5120107) q[2];
sx q[2];
rz(-1.0967983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67657802) q[1];
sx q[1];
rz(-1.4099219) q[1];
sx q[1];
rz(-1.7382998) q[1];
rz(-pi) q[2];
rz(-1.9188493) q[3];
sx q[3];
rz(-2.9677941) q[3];
sx q[3];
rz(-2.9617975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6685247) q[2];
sx q[2];
rz(-2.1418085) q[2];
sx q[2];
rz(2.0295985) q[2];
rz(2.9396577) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(0.37566617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875459) q[0];
sx q[0];
rz(-2.9534979) q[0];
sx q[0];
rz(1.4659708) q[0];
rz(-1.6000043) q[1];
sx q[1];
rz(-1.8406638) q[1];
sx q[1];
rz(-2.8864554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8153471) q[0];
sx q[0];
rz(-1.4370942) q[0];
sx q[0];
rz(0.017450138) q[0];
rz(-2.8053984) q[2];
sx q[2];
rz(-0.80226919) q[2];
sx q[2];
rz(2.4870869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7950739) q[1];
sx q[1];
rz(-0.28920214) q[1];
sx q[1];
rz(1.8498011) q[1];
x q[2];
rz(0.76948036) q[3];
sx q[3];
rz(-1.4325241) q[3];
sx q[3];
rz(-1.6618123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61350358) q[2];
sx q[2];
rz(-1.1897503) q[2];
sx q[2];
rz(0.36821723) q[2];
rz(-0.28299371) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(0.32576573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014864347) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(2.6425843) q[1];
sx q[1];
rz(-1.9022763) q[1];
sx q[1];
rz(0.3872445) q[1];
rz(-0.3327315) q[2];
sx q[2];
rz(-1.6920148) q[2];
sx q[2];
rz(0.68670338) q[2];
rz(-3.079703) q[3];
sx q[3];
rz(-1.7714942) q[3];
sx q[3];
rz(-1.9006806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
