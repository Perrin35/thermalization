OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(-0.48506919) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(5.8689868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.480455) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(-0.35868355) q[0];
x q[1];
rz(0.094751058) q[2];
sx q[2];
rz(-2.13846) q[2];
sx q[2];
rz(1.3324347) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.58798446) q[1];
sx q[1];
rz(-1.506664) q[1];
sx q[1];
rz(-1.9967805) q[1];
rz(2.8768086) q[3];
sx q[3];
rz(-1.4720819) q[3];
sx q[3];
rz(-0.36987723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(0.031575354) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-1.0715276) q[1];
sx q[1];
rz(2.6020715) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(-3.0898068) q[0];
x q[1];
rz(-2.6846634) q[2];
sx q[2];
rz(-0.77605844) q[2];
sx q[2];
rz(-1.4571783) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5121465) q[1];
sx q[1];
rz(-1.0735895) q[1];
sx q[1];
rz(0.21558233) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86124729) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(-1.4351821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4743621) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(-1.386863) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(0.2562491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2494333) q[0];
sx q[0];
rz(-1.8349577) q[0];
sx q[0];
rz(0.63412068) q[0];
rz(-pi) q[1];
rz(2.0492378) q[2];
sx q[2];
rz(-1.1698327) q[2];
sx q[2];
rz(-1.608633) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5445404) q[1];
sx q[1];
rz(-1.2433194) q[1];
sx q[1];
rz(0.87045963) q[1];
rz(-pi) q[2];
rz(-1.6085298) q[3];
sx q[3];
rz(-1.3940666) q[3];
sx q[3];
rz(0.43921525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(-0.57470542) q[2];
rz(1.3556708) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(1.2333966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.213585) q[0];
sx q[0];
rz(-1.7127697) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(2.4096699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703092) q[0];
sx q[0];
rz(-2.0844315) q[0];
sx q[0];
rz(1.0259823) q[0];
rz(-pi) q[1];
rz(0.71728431) q[2];
sx q[2];
rz(-1.4826164) q[2];
sx q[2];
rz(0.35536534) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4983474) q[1];
sx q[1];
rz(-2.8162662) q[1];
sx q[1];
rz(2.494032) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51283522) q[3];
sx q[3];
rz(-2.8898015) q[3];
sx q[3];
rz(-0.73392111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4449473) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(2.7643519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-0.60931283) q[1];
sx q[1];
rz(-0.79777065) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3482462) q[0];
sx q[0];
rz(-3.0214546) q[0];
sx q[0];
rz(1.5915464) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5227929) q[2];
sx q[2];
rz(-0.36703645) q[2];
sx q[2];
rz(0.4180846) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.55802) q[1];
sx q[1];
rz(-0.94545525) q[1];
sx q[1];
rz(-3.0261092) q[1];
rz(-pi) q[2];
rz(1.2410774) q[3];
sx q[3];
rz(-1.9922678) q[3];
sx q[3];
rz(-0.7082522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.34565869) q[2];
sx q[2];
rz(-0.63085932) q[2];
sx q[2];
rz(2.8395555) q[2];
rz(-1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-0.091826037) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(-1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(3.0715122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1228186) q[0];
sx q[0];
rz(-0.2304603) q[0];
sx q[0];
rz(-0.83968681) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3577865) q[2];
sx q[2];
rz(-1.1597826) q[2];
sx q[2];
rz(-2.59936) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97092123) q[1];
sx q[1];
rz(-1.0075924) q[1];
sx q[1];
rz(1.7621653) q[1];
x q[2];
rz(3.1061884) q[3];
sx q[3];
rz(-1.4962713) q[3];
sx q[3];
rz(-0.41985928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30248102) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(-2.5682209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(-2.281718) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-0.87221968) q[1];
sx q[1];
rz(3.133657) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82848362) q[0];
sx q[0];
rz(-2.5223753) q[0];
sx q[0];
rz(2.6909268) q[0];
rz(0.43086149) q[2];
sx q[2];
rz(-0.92509809) q[2];
sx q[2];
rz(2.8889887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.24212813) q[1];
sx q[1];
rz(-2.2367396) q[1];
sx q[1];
rz(2.4813586) q[1];
x q[2];
rz(-0.90585917) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(2.0743899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.8442644) q[3];
sx q[3];
rz(0.90464512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96173441) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-2.7767048) q[0];
rz(0.94003135) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(-1.6392802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84425981) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(0.25015932) q[0];
rz(-pi) q[1];
rz(-2.6697568) q[2];
sx q[2];
rz(-0.89315692) q[2];
sx q[2];
rz(-3.133528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9099721) q[1];
sx q[1];
rz(-2.2044704) q[1];
sx q[1];
rz(0.85068591) q[1];
rz(1.1542529) q[3];
sx q[3];
rz(-2.917218) q[3];
sx q[3];
rz(-0.32240543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6442948) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(-1.0650744) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97312462) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(0.43564963) q[0];
rz(1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(-0.41697821) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1736261) q[0];
sx q[0];
rz(-0.91714232) q[0];
sx q[0];
rz(-0.047441479) q[0];
rz(0.66395335) q[2];
sx q[2];
rz(-1.10154) q[2];
sx q[2];
rz(-2.9479153) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.13874395) q[1];
sx q[1];
rz(-1.7094304) q[1];
sx q[1];
rz(1.0671875) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3973893) q[3];
sx q[3];
rz(-0.19072285) q[3];
sx q[3];
rz(0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(-0.22582516) q[2];
rz(2.9337692) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(0.58661714) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.5243994) q[0];
rz(2.1879451) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.8189925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029862558) q[0];
sx q[0];
rz(-2.7435281) q[0];
sx q[0];
rz(-1.2533305) q[0];
rz(-1.0820504) q[2];
sx q[2];
rz(-0.49968038) q[2];
sx q[2];
rz(2.3479455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0304071) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(-2.2305957) q[1];
rz(-pi) q[2];
rz(-1.3721458) q[3];
sx q[3];
rz(-2.47654) q[3];
sx q[3];
rz(-2.8952451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(1.0160758) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(-0.55541742) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(-1.7779508) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(0.017756391) q[2];
sx q[2];
rz(-2.6323071) q[2];
sx q[2];
rz(1.4321362) q[2];
rz(-0.090311269) q[3];
sx q[3];
rz(-1.8009381) q[3];
sx q[3];
rz(1.2931852) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];