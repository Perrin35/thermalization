OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(-0.68343502) q[0];
sx q[0];
rz(0.47877065) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(-0.64136139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.1874275) q[0];
sx q[0];
rz(-1.8281561) q[0];
rz(-pi) q[1];
rz(2.7641069) q[2];
sx q[2];
rz(-1.1698206) q[2];
sx q[2];
rz(0.12667835) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76170834) q[1];
sx q[1];
rz(-2.0239081) q[1];
sx q[1];
rz(2.3945827) q[1];
rz(-pi) q[2];
rz(-0.62946837) q[3];
sx q[3];
rz(-1.9068204) q[3];
sx q[3];
rz(-2.0365086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59149867) q[2];
sx q[2];
rz(-1.9248362) q[2];
sx q[2];
rz(-0.47810289) q[2];
rz(-1.6889307) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(1.0189198) q[0];
rz(-1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(2.5085124) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7019254) q[0];
sx q[0];
rz(-0.81322008) q[0];
sx q[0];
rz(1.4600091) q[0];
rz(-pi) q[1];
rz(-2.1585629) q[2];
sx q[2];
rz(-2.2289742) q[2];
sx q[2];
rz(-2.7764729) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9620348) q[1];
sx q[1];
rz(-1.9704559) q[1];
sx q[1];
rz(2.7707997) q[1];
rz(-0.31140621) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(2.1962375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35976609) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(-0.11745545) q[2];
rz(-2.84058) q[3];
sx q[3];
rz(-1.8071226) q[3];
sx q[3];
rz(1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(-0.088407956) q[0];
rz(2.03405) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-0.69082469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3120964) q[0];
sx q[0];
rz(-1.4179686) q[0];
sx q[0];
rz(1.6468847) q[0];
rz(-pi) q[1];
rz(1.3356528) q[2];
sx q[2];
rz(-2.3217839) q[2];
sx q[2];
rz(0.90607925) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.01268) q[1];
sx q[1];
rz(-1.3887654) q[1];
sx q[1];
rz(-3.1151031) q[1];
x q[2];
rz(2.7615943) q[3];
sx q[3];
rz(-1.899154) q[3];
sx q[3];
rz(-2.2896555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.92418015) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(-0.10647354) q[2];
rz(-1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(0.52247125) q[0];
rz(-2.8126295) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(0.55363384) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75008167) q[0];
sx q[0];
rz(-0.9010074) q[0];
sx q[0];
rz(-2.1704587) q[0];
rz(-pi) q[1];
x q[1];
rz(2.373898) q[2];
sx q[2];
rz(-0.97937095) q[2];
sx q[2];
rz(-2.5007345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1045038) q[1];
sx q[1];
rz(-1.4375086) q[1];
sx q[1];
rz(-0.78661211) q[1];
rz(-pi) q[2];
rz(-0.028285154) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(2.459211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(2.6468357) q[2];
rz(-0.90302145) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(-1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49098) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(2.0565128) q[0];
rz(-0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(-0.18403149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9404011) q[0];
sx q[0];
rz(-1.860025) q[0];
sx q[0];
rz(-1.679923) q[0];
rz(-0.63921914) q[2];
sx q[2];
rz(-2.6972065) q[2];
sx q[2];
rz(-1.5612615) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6726471) q[1];
sx q[1];
rz(-2.1429859) q[1];
sx q[1];
rz(-2.1395626) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3695413) q[3];
sx q[3];
rz(-2.2866837) q[3];
sx q[3];
rz(1.2224689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.90157834) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(1.1408172) q[2];
rz(-0.42282894) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(2.254803) q[0];
rz(-2.9011762) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(2.9930847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7667023) q[0];
sx q[0];
rz(-2.3935211) q[0];
sx q[0];
rz(0.23776777) q[0];
rz(-pi) q[1];
rz(-2.912942) q[2];
sx q[2];
rz(-1.4637404) q[2];
sx q[2];
rz(0.46496898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0420694) q[1];
sx q[1];
rz(-1.6798786) q[1];
sx q[1];
rz(-1.8314349) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6753747) q[3];
sx q[3];
rz(-1.2122452) q[3];
sx q[3];
rz(-2.8311604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85577661) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(-2.6532069) q[2];
rz(2.6521902) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(-2.3235902) q[0];
rz(1.2524293) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-1.12524) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.274652) q[0];
sx q[0];
rz(-0.38892239) q[0];
sx q[0];
rz(1.2961943) q[0];
rz(-pi) q[1];
rz(-1.2721887) q[2];
sx q[2];
rz(-1.719279) q[2];
sx q[2];
rz(-2.6872925) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.352467) q[1];
sx q[1];
rz(-1.2965634) q[1];
sx q[1];
rz(1.2624361) q[1];
rz(-pi) q[2];
rz(-2.1358228) q[3];
sx q[3];
rz(-2.1566448) q[3];
sx q[3];
rz(2.1573531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(-1.6623496) q[3];
sx q[3];
rz(-2.1846266) q[3];
sx q[3];
rz(-1.7361599) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97061625) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(1.6059426) q[0];
rz(-1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-1.01064) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7692524) q[0];
sx q[0];
rz(-1.7307889) q[0];
sx q[0];
rz(0.82323797) q[0];
x q[1];
rz(2.8433617) q[2];
sx q[2];
rz(-1.094162) q[2];
sx q[2];
rz(-0.44520608) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6992221) q[1];
sx q[1];
rz(-1.7915627) q[1];
sx q[1];
rz(-1.8063596) q[1];
rz(-pi) q[2];
rz(0.2145433) q[3];
sx q[3];
rz(-0.40896591) q[3];
sx q[3];
rz(-1.6681125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6179787) q[2];
sx q[2];
rz(-0.31948677) q[2];
sx q[2];
rz(0.69407216) q[2];
rz(-0.56898919) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(-2.2038961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-3.0555608) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(-0.49474299) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(3.0659952) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6583017) q[0];
sx q[0];
rz(-1.0169944) q[0];
sx q[0];
rz(2.7730586) q[0];
rz(-pi) q[1];
rz(1.2443301) q[2];
sx q[2];
rz(-0.18507659) q[2];
sx q[2];
rz(-1.2125804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96502111) q[1];
sx q[1];
rz(-0.88639835) q[1];
sx q[1];
rz(-1.4375163) q[1];
rz(-pi) q[2];
rz(0.62065403) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.8295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8081234) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(-1.3405651) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(0.17679581) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(-0.44100824) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1045751) q[0];
sx q[0];
rz(-1.772176) q[0];
sx q[0];
rz(1.7778394) q[0];
rz(-pi) q[1];
x q[1];
rz(2.987864) q[2];
sx q[2];
rz(-1.4970461) q[2];
sx q[2];
rz(-1.7388625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1383789) q[1];
sx q[1];
rz(-0.56204501) q[1];
sx q[1];
rz(-0.89488645) q[1];
rz(-1.2714083) q[3];
sx q[3];
rz(-2.1736439) q[3];
sx q[3];
rz(0.30764461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(2.8397172) q[2];
rz(2.2484696) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72398913) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(-0.026731116) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-0.46840657) q[2];
sx q[2];
rz(-0.9006587) q[2];
sx q[2];
rz(-2.547154) q[2];
rz(-1.0711014) q[3];
sx q[3];
rz(-1.0715967) q[3];
sx q[3];
rz(-1.788492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
