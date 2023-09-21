OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5972714) q[0];
sx q[0];
rz(-0.40364021) q[0];
sx q[0];
rz(-0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(-1.4189781) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0231853) q[0];
sx q[0];
rz(-1.7405131) q[0];
sx q[0];
rz(-2.8939308) q[0];
rz(-1.8364036) q[2];
sx q[2];
rz(-2.0678389) q[2];
sx q[2];
rz(-3.1211987) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7367243) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(3.0800372) q[1];
rz(-pi) q[2];
rz(1.2526413) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(0.51154256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(0.35787004) q[2];
rz(-0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(0.068280846) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(-0.65111792) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40819528) q[0];
sx q[0];
rz(-1.3696284) q[0];
sx q[0];
rz(0.64543076) q[0];
rz(-pi) q[1];
rz(0.5257734) q[2];
sx q[2];
rz(-1.5915046) q[2];
sx q[2];
rz(-2.6244342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5755641) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(-2.6355987) q[1];
rz(-2.7495456) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(0.80313659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-0.90332705) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(-2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(2.6480411) q[0];
rz(1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-1.0292056) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632336) q[0];
sx q[0];
rz(-1.2326816) q[0];
sx q[0];
rz(2.6792206) q[0];
x q[1];
rz(-0.87666338) q[2];
sx q[2];
rz(-2.6138517) q[2];
sx q[2];
rz(-1.1610247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.3612559) q[1];
sx q[1];
rz(0.31063147) q[1];
x q[2];
rz(-2.1152705) q[3];
sx q[3];
rz(-1.4635651) q[3];
sx q[3];
rz(0.925392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-2.4659618) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-0.40859616) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-2.2600007) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170167) q[0];
sx q[0];
rz(-1.3991742) q[0];
sx q[0];
rz(0.73715985) q[0];
rz(-0.84882952) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(-1.3918849) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2320764) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(2.3823649) q[1];
rz(-pi) q[2];
rz(-1.5041755) q[3];
sx q[3];
rz(-1.1542218) q[3];
sx q[3];
rz(1.6601603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9920138) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-2.8240906) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(0.11725765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7754585) q[0];
sx q[0];
rz(-1.701208) q[0];
sx q[0];
rz(0.58701697) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3567032) q[2];
sx q[2];
rz(-0.97852409) q[2];
sx q[2];
rz(-1.585373) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.458141) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(1.7276006) q[1];
rz(0.076905964) q[3];
sx q[3];
rz(-2.1566026) q[3];
sx q[3];
rz(2.2211071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7964898) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538486) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(3.0336753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6861434) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(2.7129052) q[0];
rz(-pi) q[1];
rz(-3.1015322) q[2];
sx q[2];
rz(-1.5594348) q[2];
sx q[2];
rz(-2.5324015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94090998) q[1];
sx q[1];
rz(-0.74146491) q[1];
sx q[1];
rz(-2.7156668) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88166724) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(-1.4985639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-1.1090013) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(-3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7727707) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(-2.7203454) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(-0.79963911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0945471) q[0];
sx q[0];
rz(-0.73091113) q[0];
sx q[0];
rz(2.7546309) q[0];
rz(1.854935) q[2];
sx q[2];
rz(-1.8737027) q[2];
sx q[2];
rz(-2.6877833) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0529424) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(0.36693962) q[1];
rz(-pi) q[2];
rz(-2.1361254) q[3];
sx q[3];
rz(-1.4544832) q[3];
sx q[3];
rz(-3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(2.183765) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9668982) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(2.877537) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(-2.82428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89462751) q[0];
sx q[0];
rz(-1.4980226) q[0];
sx q[0];
rz(-1.9586246) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9192341) q[2];
sx q[2];
rz(-1.7218044) q[2];
sx q[2];
rz(0.82049557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4006153) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(0.11465794) q[1];
rz(-pi) q[2];
rz(-1.4044912) q[3];
sx q[3];
rz(-2.2095223) q[3];
sx q[3];
rz(-0.8730264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43549609) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(1.502011) q[2];
rz(-2.8912985) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.1423473) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89649993) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(1.7171575) q[0];
rz(2.6622488) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(0.11553484) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53287017) q[0];
sx q[0];
rz(-2.1358129) q[0];
sx q[0];
rz(0.67824058) q[0];
rz(-pi) q[1];
rz(-0.71474448) q[2];
sx q[2];
rz(-1.359316) q[2];
sx q[2];
rz(-1.5541058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9563453) q[1];
sx q[1];
rz(-0.51925175) q[1];
sx q[1];
rz(-0.96038702) q[1];
x q[2];
rz(-1.450591) q[3];
sx q[3];
rz(-0.87814858) q[3];
sx q[3];
rz(1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7675161) q[0];
sx q[0];
rz(-2.9449468) q[0];
sx q[0];
rz(-1.0023414) q[0];
rz(-pi) q[1];
rz(-1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(-1.64738) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1684432) q[1];
sx q[1];
rz(-1.397555) q[1];
sx q[1];
rz(-3.1107305) q[1];
x q[2];
rz(1.2606603) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(-0.96084259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.6258378) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(2.1102171) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2172858) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-2.968593) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(1.8299673) q[3];
sx q[3];
rz(-1.6298686) q[3];
sx q[3];
rz(2.009404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
