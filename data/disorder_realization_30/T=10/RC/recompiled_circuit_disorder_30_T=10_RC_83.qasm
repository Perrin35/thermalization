OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6501939) q[0];
sx q[0];
rz(-2.8770652) q[0];
sx q[0];
rz(-2.7471623) q[0];
rz(0.0061622942) q[1];
sx q[1];
rz(-0.34024629) q[1];
sx q[1];
rz(1.9415829) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7317176) q[0];
sx q[0];
rz(-3.0933676) q[0];
sx q[0];
rz(1.6884786) q[0];
rz(0.94304396) q[2];
sx q[2];
rz(-2.606751) q[2];
sx q[2];
rz(-0.56378555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8438699) q[1];
sx q[1];
rz(-0.6134609) q[1];
sx q[1];
rz(-2.2754199) q[1];
rz(0.37813152) q[3];
sx q[3];
rz(-1.5298801) q[3];
sx q[3];
rz(3.0331628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3068984) q[2];
sx q[2];
rz(-1.4935741) q[2];
sx q[2];
rz(-2.8519894) q[2];
rz(-0.87537193) q[3];
sx q[3];
rz(-2.1353728) q[3];
sx q[3];
rz(-3.0818821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.61525476) q[0];
sx q[0];
rz(-0.86741388) q[0];
sx q[0];
rz(-0.34399024) q[0];
rz(3.0572609) q[1];
sx q[1];
rz(-2.4721959) q[1];
sx q[1];
rz(1.3551691) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7536613) q[0];
sx q[0];
rz(-0.61457115) q[0];
sx q[0];
rz(1.4580926) q[0];
rz(-pi) q[1];
rz(-0.64135735) q[2];
sx q[2];
rz(-2.9396082) q[2];
sx q[2];
rz(0.55822492) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.720571) q[1];
sx q[1];
rz(-2.6877626) q[1];
sx q[1];
rz(-1.4245016) q[1];
rz(-pi) q[2];
rz(3.063952) q[3];
sx q[3];
rz(-1.3878126) q[3];
sx q[3];
rz(1.3444572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29558674) q[2];
sx q[2];
rz(-2.3082374) q[2];
sx q[2];
rz(-2.6039092) q[2];
rz(-0.50283557) q[3];
sx q[3];
rz(-2.7167795) q[3];
sx q[3];
rz(-1.0104377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780592) q[0];
sx q[0];
rz(-2.6038267) q[0];
sx q[0];
rz(0.64087254) q[0];
rz(2.3928941) q[1];
sx q[1];
rz(-1.0083895) q[1];
sx q[1];
rz(2.0764988) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070806064) q[0];
sx q[0];
rz(-2.7822128) q[0];
sx q[0];
rz(-1.552836) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58653411) q[2];
sx q[2];
rz(-2.1792992) q[2];
sx q[2];
rz(1.3003132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0712191) q[1];
sx q[1];
rz(-0.95446903) q[1];
sx q[1];
rz(3.0973869) q[1];
rz(2.2933526) q[3];
sx q[3];
rz(-1.3812997) q[3];
sx q[3];
rz(3.0998067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.76434) q[2];
sx q[2];
rz(-1.9621907) q[2];
sx q[2];
rz(-2.4242145) q[2];
rz(2.453089) q[3];
sx q[3];
rz(-0.62763667) q[3];
sx q[3];
rz(-2.8440516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82729572) q[0];
sx q[0];
rz(-1.2001487) q[0];
sx q[0];
rz(-2.8919019) q[0];
rz(1.0149792) q[1];
sx q[1];
rz(-2.8453638) q[1];
sx q[1];
rz(0.011118523) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5549094) q[0];
sx q[0];
rz(-1.1197829) q[0];
sx q[0];
rz(0.45278544) q[0];
x q[1];
rz(-0.87984933) q[2];
sx q[2];
rz(-0.86266154) q[2];
sx q[2];
rz(-0.81530064) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1383458) q[1];
sx q[1];
rz(-1.0838638) q[1];
sx q[1];
rz(2.7702615) q[1];
x q[2];
rz(0.31828493) q[3];
sx q[3];
rz(-1.3878229) q[3];
sx q[3];
rz(1.1312248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49542385) q[2];
sx q[2];
rz(-1.255722) q[2];
sx q[2];
rz(0.28309506) q[2];
rz(2.4781573) q[3];
sx q[3];
rz(-2.5484271) q[3];
sx q[3];
rz(-2.2535113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11113142) q[0];
sx q[0];
rz(-2.9009394) q[0];
sx q[0];
rz(-2.8097613) q[0];
rz(-0.49452531) q[1];
sx q[1];
rz(-1.8437513) q[1];
sx q[1];
rz(1.3269075) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622612) q[0];
sx q[0];
rz(-2.3752897) q[0];
sx q[0];
rz(-2.9192231) q[0];
rz(-pi) q[1];
rz(-0.58381501) q[2];
sx q[2];
rz(-2.2540255) q[2];
sx q[2];
rz(0.77606397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9304282) q[1];
sx q[1];
rz(-1.8976364) q[1];
sx q[1];
rz(1.8153166) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96552403) q[3];
sx q[3];
rz(-11*pi/13) q[3];
sx q[3];
rz(2.0616639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.999324) q[2];
sx q[2];
rz(-0.48196718) q[2];
sx q[2];
rz(-1.2456606) q[2];
rz(1.8866395) q[3];
sx q[3];
rz(-2.3001223) q[3];
sx q[3];
rz(0.83827034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6923043) q[0];
sx q[0];
rz(-0.0031539991) q[0];
sx q[0];
rz(-0.6814878) q[0];
rz(-0.20755126) q[1];
sx q[1];
rz(-0.46336585) q[1];
sx q[1];
rz(2.025827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078243144) q[0];
sx q[0];
rz(-2.3176498) q[0];
sx q[0];
rz(1.0043762) q[0];
rz(-pi) q[1];
rz(2.0954779) q[2];
sx q[2];
rz(-2.4002053) q[2];
sx q[2];
rz(1.807883) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9876154) q[1];
sx q[1];
rz(-2.1591641) q[1];
sx q[1];
rz(-1.5495367) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3917771) q[3];
sx q[3];
rz(-0.90913032) q[3];
sx q[3];
rz(-0.20565198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21268022) q[2];
sx q[2];
rz(-1.2968411) q[2];
sx q[2];
rz(0.35432717) q[2];
rz(-2.8220693) q[3];
sx q[3];
rz(-1.9890246) q[3];
sx q[3];
rz(0.37187809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324683) q[0];
sx q[0];
rz(-0.23129825) q[0];
sx q[0];
rz(-0.67434597) q[0];
rz(2.0293503) q[1];
sx q[1];
rz(-2.4770885) q[1];
sx q[1];
rz(-0.56232125) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494941) q[0];
sx q[0];
rz(-1.5559762) q[0];
sx q[0];
rz(-0.90256079) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33195915) q[2];
sx q[2];
rz(-1.6662285) q[2];
sx q[2];
rz(-1.7711668) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6866236) q[1];
sx q[1];
rz(-1.3775871) q[1];
sx q[1];
rz(-1.4229694) q[1];
x q[2];
rz(1.698877) q[3];
sx q[3];
rz(-2.5411385) q[3];
sx q[3];
rz(-1.2680935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0397296) q[2];
sx q[2];
rz(-0.9938643) q[2];
sx q[2];
rz(2.8015461) q[2];
rz(0.21751054) q[3];
sx q[3];
rz(-0.9483996) q[3];
sx q[3];
rz(-0.31869179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44889221) q[0];
sx q[0];
rz(-0.046148766) q[0];
sx q[0];
rz(-0.39644077) q[0];
rz(3.0026644) q[1];
sx q[1];
rz(-0.46008343) q[1];
sx q[1];
rz(1.6202392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053145807) q[0];
sx q[0];
rz(-1.7674812) q[0];
sx q[0];
rz(1.501207) q[0];
x q[1];
rz(-2.7395758) q[2];
sx q[2];
rz(-2.0094299) q[2];
sx q[2];
rz(-1.5356262) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.13008598) q[1];
sx q[1];
rz(-2.3475921) q[1];
sx q[1];
rz(1.4873051) q[1];
x q[2];
rz(1.7117386) q[3];
sx q[3];
rz(-2.4574276) q[3];
sx q[3];
rz(1.14389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56269318) q[2];
sx q[2];
rz(-2.0843299) q[2];
sx q[2];
rz(-2.8472624) q[2];
rz(-2.0108022) q[3];
sx q[3];
rz(-1.7629938) q[3];
sx q[3];
rz(0.99564266) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49333736) q[0];
sx q[0];
rz(-0.91518891) q[0];
sx q[0];
rz(-0.35933581) q[0];
rz(-0.94611478) q[1];
sx q[1];
rz(-2.7455536) q[1];
sx q[1];
rz(2.8709581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7662738) q[0];
sx q[0];
rz(-3.0941071) q[0];
sx q[0];
rz(2.4256698) q[0];
rz(-pi) q[1];
rz(1.0893906) q[2];
sx q[2];
rz(-0.82197661) q[2];
sx q[2];
rz(-0.35441986) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7120413) q[1];
sx q[1];
rz(-1.5263285) q[1];
sx q[1];
rz(0.045750381) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1569096) q[3];
sx q[3];
rz(-2.1911748) q[3];
sx q[3];
rz(-1.2445205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82751194) q[2];
sx q[2];
rz(-2.9120047) q[2];
sx q[2];
rz(2.6861526) q[2];
rz(0.81196249) q[3];
sx q[3];
rz(-2.2596695) q[3];
sx q[3];
rz(0.5493831) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51046002) q[0];
sx q[0];
rz(-1.5083418) q[0];
sx q[0];
rz(0.73927885) q[0];
rz(2.9108858) q[1];
sx q[1];
rz(-0.46805996) q[1];
sx q[1];
rz(-0.49490067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65942818) q[0];
sx q[0];
rz(-2.6539301) q[0];
sx q[0];
rz(1.6517261) q[0];
rz(0.020936326) q[2];
sx q[2];
rz(-2.8136721) q[2];
sx q[2];
rz(-2.7431969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7086664) q[1];
sx q[1];
rz(-0.64879829) q[1];
sx q[1];
rz(1.385958) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4570191) q[3];
sx q[3];
rz(-2.4416231) q[3];
sx q[3];
rz(1.1752807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0080002) q[2];
sx q[2];
rz(-2.0501037) q[2];
sx q[2];
rz(0.32785329) q[2];
rz(-0.19206364) q[3];
sx q[3];
rz(-0.24774194) q[3];
sx q[3];
rz(2.1081934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(1.0192169) q[0];
sx q[0];
rz(-1.6288971) q[0];
sx q[0];
rz(-1.6678641) q[0];
rz(2.3090251) q[1];
sx q[1];
rz(-1.4214129) q[1];
sx q[1];
rz(1.3690154) q[1];
rz(-2.9359948) q[2];
sx q[2];
rz(-1.0369221) q[2];
sx q[2];
rz(2.592929) q[2];
rz(1.8264063) q[3];
sx q[3];
rz(-1.3127463) q[3];
sx q[3];
rz(-1.128872) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
