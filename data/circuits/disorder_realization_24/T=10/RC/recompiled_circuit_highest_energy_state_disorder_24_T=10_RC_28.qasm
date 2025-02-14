OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10676323) q[0];
sx q[0];
rz(-1.8101298) q[0];
sx q[0];
rz(-1.3504299) q[0];
rz(-2.8332233) q[1];
sx q[1];
rz(-2.8928533) q[1];
sx q[1];
rz(-2.1623936) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0686183) q[0];
sx q[0];
rz(-0.81356293) q[0];
sx q[0];
rz(-1.6723775) q[0];
x q[1];
rz(-2.2297284) q[2];
sx q[2];
rz(-1.2871337) q[2];
sx q[2];
rz(2.5549169) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56899988) q[1];
sx q[1];
rz(-1.6385983) q[1];
sx q[1];
rz(2.7876233) q[1];
x q[2];
rz(1.3967674) q[3];
sx q[3];
rz(-1.9367894) q[3];
sx q[3];
rz(-0.39453615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(2.5085874) q[2];
rz(-0.64166075) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(-2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68384701) q[0];
sx q[0];
rz(-2.1774543) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(-1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(0.32018426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31486785) q[0];
sx q[0];
rz(-1.9786252) q[0];
sx q[0];
rz(2.639222) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9302017) q[2];
sx q[2];
rz(-2.0703982) q[2];
sx q[2];
rz(1.4818986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.089848891) q[1];
sx q[1];
rz(-1.4898572) q[1];
sx q[1];
rz(-1.5481987) q[1];
rz(-pi) q[2];
rz(-3.0820936) q[3];
sx q[3];
rz(-1.3493285) q[3];
sx q[3];
rz(2.6513316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1227293) q[2];
sx q[2];
rz(-1.2168695) q[2];
sx q[2];
rz(-1.8625721) q[2];
rz(-0.095712885) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(-1.3191282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1529481) q[0];
sx q[0];
rz(-2.4929292) q[0];
sx q[0];
rz(-2.1107819) q[0];
rz(-0.55054322) q[1];
sx q[1];
rz(-0.85854733) q[1];
sx q[1];
rz(1.8969089) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36280945) q[0];
sx q[0];
rz(-0.16805695) q[0];
sx q[0];
rz(-0.3248949) q[0];
x q[1];
rz(0.16764201) q[2];
sx q[2];
rz(-0.9603921) q[2];
sx q[2];
rz(-0.33627015) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.74725973) q[1];
sx q[1];
rz(-1.7752517) q[1];
sx q[1];
rz(2.2482199) q[1];
rz(1.9774417) q[3];
sx q[3];
rz(-2.3465524) q[3];
sx q[3];
rz(0.94889489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9024258) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(1.172056) q[2];
rz(-0.82550448) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(2.4765305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.091752) q[0];
sx q[0];
rz(-1.4565775) q[0];
sx q[0];
rz(-2.7276584) q[0];
rz(-2.6992758) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(0.54534674) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6028311) q[0];
sx q[0];
rz(-1.6206546) q[0];
sx q[0];
rz(0.2754579) q[0];
rz(2.2121934) q[2];
sx q[2];
rz(-2.2365733) q[2];
sx q[2];
rz(-2.5432472) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2703823) q[1];
sx q[1];
rz(-0.59894136) q[1];
sx q[1];
rz(-1.0652871) q[1];
rz(-pi) q[2];
rz(0.14800565) q[3];
sx q[3];
rz(-1.0665575) q[3];
sx q[3];
rz(-1.3324236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44117323) q[2];
sx q[2];
rz(-1.8808695) q[2];
sx q[2];
rz(2.9803661) q[2];
rz(1.5876596) q[3];
sx q[3];
rz(-2.5416538) q[3];
sx q[3];
rz(0.59066311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081838354) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(0.72342122) q[0];
rz(-1.3149892) q[1];
sx q[1];
rz(-1.864121) q[1];
sx q[1];
rz(-2.1037197) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53053601) q[0];
sx q[0];
rz(-0.7270455) q[0];
sx q[0];
rz(-0.1188936) q[0];
rz(-2.5332301) q[2];
sx q[2];
rz(-0.82348292) q[2];
sx q[2];
rz(-1.057511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65601634) q[1];
sx q[1];
rz(-0.52158251) q[1];
sx q[1];
rz(-2.5822469) q[1];
x q[2];
rz(2.5955022) q[3];
sx q[3];
rz(-1.4442863) q[3];
sx q[3];
rz(1.5777335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3967241) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(-0.57207668) q[2];
rz(2.5165596) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(-1.2264138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142757) q[0];
sx q[0];
rz(-0.033997424) q[0];
sx q[0];
rz(-2.6183364) q[0];
rz(-1.8190067) q[1];
sx q[1];
rz(-1.4679642) q[1];
sx q[1];
rz(2.5214213) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70411982) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(1.062514) q[0];
x q[1];
rz(0.1441346) q[2];
sx q[2];
rz(-0.47283334) q[2];
sx q[2];
rz(-2.0236058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6002308) q[1];
sx q[1];
rz(-2.5087452) q[1];
sx q[1];
rz(-1.71943) q[1];
rz(2.9240588) q[3];
sx q[3];
rz(-1.7420046) q[3];
sx q[3];
rz(-1.8678766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72914499) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(1.8577417) q[2];
rz(0.9681975) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(1.5549972) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17484434) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(-1.211776) q[0];
rz(-0.9696331) q[1];
sx q[1];
rz(-1.8200834) q[1];
sx q[1];
rz(-1.2744354) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6912708) q[0];
sx q[0];
rz(-0.88435882) q[0];
sx q[0];
rz(-2.9458988) q[0];
rz(-1.6658804) q[2];
sx q[2];
rz(-2.3543752) q[2];
sx q[2];
rz(-1.9958391) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.695622) q[1];
sx q[1];
rz(-1.8730764) q[1];
sx q[1];
rz(0.92960371) q[1];
rz(0.46044965) q[3];
sx q[3];
rz(-1.2795953) q[3];
sx q[3];
rz(1.2429383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6756639) q[2];
sx q[2];
rz(-1.6589386) q[2];
sx q[2];
rz(-0.3234123) q[2];
rz(0.0023500738) q[3];
sx q[3];
rz(-1.4582783) q[3];
sx q[3];
rz(0.6215483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6482786) q[0];
sx q[0];
rz(-1.6596153) q[0];
sx q[0];
rz(1.3508654) q[0];
rz(1.3510652) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(1.2877119) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0553594) q[0];
sx q[0];
rz(-1.6882734) q[0];
sx q[0];
rz(1.286195) q[0];
rz(-pi) q[1];
rz(1.9016198) q[2];
sx q[2];
rz(-0.74218732) q[2];
sx q[2];
rz(-2.0798707) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9629891) q[1];
sx q[1];
rz(-1.7409391) q[1];
sx q[1];
rz(-0.32802204) q[1];
rz(-pi) q[2];
rz(1.5896371) q[3];
sx q[3];
rz(-1.9196438) q[3];
sx q[3];
rz(-0.89295372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6172341) q[2];
sx q[2];
rz(-1.989216) q[2];
sx q[2];
rz(0.99346811) q[2];
rz(-1.0348381) q[3];
sx q[3];
rz(-0.60641685) q[3];
sx q[3];
rz(-0.16916999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.8096823) q[0];
sx q[0];
rz(-1.1746291) q[0];
sx q[0];
rz(-1.6140953) q[0];
rz(0.88036674) q[1];
sx q[1];
rz(-2.7208734) q[1];
sx q[1];
rz(-2.8299433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0880222) q[0];
sx q[0];
rz(-0.35614518) q[0];
sx q[0];
rz(-0.97266622) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3394322) q[2];
sx q[2];
rz(-1.5867434) q[2];
sx q[2];
rz(-2.0248264) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2276409) q[1];
sx q[1];
rz(-2.3415509) q[1];
sx q[1];
rz(-2.6197919) q[1];
x q[2];
rz(1.9945456) q[3];
sx q[3];
rz(-0.69931385) q[3];
sx q[3];
rz(2.5831403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46939048) q[2];
sx q[2];
rz(-0.92350525) q[2];
sx q[2];
rz(-1.3908609) q[2];
rz(1.6145128) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(-1.0745777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-2.5020318) q[0];
sx q[0];
rz(-1.1447516) q[0];
sx q[0];
rz(0.61258739) q[0];
rz(-2.1514429) q[1];
sx q[1];
rz(-1.1724816) q[1];
sx q[1];
rz(1.6506857) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6263437) q[0];
sx q[0];
rz(-0.83413163) q[0];
sx q[0];
rz(2.5596928) q[0];
x q[1];
rz(-0.91103001) q[2];
sx q[2];
rz(-1.3675895) q[2];
sx q[2];
rz(3.0958297) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4348169) q[1];
sx q[1];
rz(-1.5770327) q[1];
sx q[1];
rz(1.4093136) q[1];
x q[2];
rz(-0.29514031) q[3];
sx q[3];
rz(-2.4217229) q[3];
sx q[3];
rz(0.9269549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7154197) q[2];
sx q[2];
rz(-2.3121068) q[2];
sx q[2];
rz(-3.1066011) q[2];
rz(-1.70111) q[3];
sx q[3];
rz(-2.248843) q[3];
sx q[3];
rz(2.6006202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69915019) q[0];
sx q[0];
rz(-0.89542605) q[0];
sx q[0];
rz(-0.17740346) q[0];
rz(-1.9433446) q[1];
sx q[1];
rz(-1.5180963) q[1];
sx q[1];
rz(0.95388283) q[1];
rz(2.4631029) q[2];
sx q[2];
rz(-1.0889006) q[2];
sx q[2];
rz(2.0550645) q[2];
rz(-2.0318326) q[3];
sx q[3];
rz(-2.2363935) q[3];
sx q[3];
rz(1.4105894) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
