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
rz(-2.1751997) q[0];
sx q[0];
rz(-1.3858495) q[0];
sx q[0];
rz(0.59544271) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(2.5505677) q[1];
sx q[1];
rz(12.443065) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6791413) q[0];
sx q[0];
rz(-1.9872905) q[0];
sx q[0];
rz(-0.97521675) q[0];
rz(-pi) q[1];
rz(2.3507471) q[2];
sx q[2];
rz(-2.186785) q[2];
sx q[2];
rz(0.14033422) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6001624) q[1];
sx q[1];
rz(-1.7847536) q[1];
sx q[1];
rz(0.096264953) q[1];
rz(0.38148613) q[3];
sx q[3];
rz(-0.39159989) q[3];
sx q[3];
rz(2.7158383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3143602) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(0.11191351) q[2];
rz(1.7498451) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(2.8351423) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7258485) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(2.9916905) q[0];
rz(-0.28597486) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(1.1289977) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647346) q[0];
sx q[0];
rz(-1.5405415) q[0];
sx q[0];
rz(-1.8591575) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3311884) q[2];
sx q[2];
rz(-3.0802266) q[2];
sx q[2];
rz(0.97739109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4399229) q[1];
sx q[1];
rz(-1.4050802) q[1];
sx q[1];
rz(-2.7782337) q[1];
x q[2];
rz(-2.1993447) q[3];
sx q[3];
rz(-2.0671004) q[3];
sx q[3];
rz(-1.7402349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67843208) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(-1.1506259) q[2];
rz(1.1848508) q[3];
sx q[3];
rz(-1.4169644) q[3];
sx q[3];
rz(3.1209893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612369) q[0];
sx q[0];
rz(-0.24599563) q[0];
sx q[0];
rz(2.7599957) q[0];
rz(1.8474139) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(2.9433184) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9851882) q[0];
sx q[0];
rz(-1.3964147) q[0];
sx q[0];
rz(-1.8013493) q[0];
rz(-pi) q[1];
rz(1.4691822) q[2];
sx q[2];
rz(-0.67199113) q[2];
sx q[2];
rz(2.4198614) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5845909) q[1];
sx q[1];
rz(-1.3583359) q[1];
sx q[1];
rz(-0.23834385) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1214633) q[3];
sx q[3];
rz(-1.745589) q[3];
sx q[3];
rz(-2.54769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3759489) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(-0.51885968) q[2];
rz(1.8910003) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4984703) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(-0.44922391) q[0];
rz(1.4462224) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(0.62072388) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45263824) q[0];
sx q[0];
rz(-1.7151378) q[0];
sx q[0];
rz(2.4792433) q[0];
rz(0.63346699) q[2];
sx q[2];
rz(-1.6704428) q[2];
sx q[2];
rz(1.5657305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6606969) q[1];
sx q[1];
rz(-0.45240739) q[1];
sx q[1];
rz(-2.6651732) q[1];
x q[2];
rz(-1.4820547) q[3];
sx q[3];
rz(-2.6013298) q[3];
sx q[3];
rz(2.1335354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(0.16109666) q[2];
rz(2.7437239) q[3];
sx q[3];
rz(-1.4784003) q[3];
sx q[3];
rz(-2.9919992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75905269) q[0];
sx q[0];
rz(-1.7790786) q[0];
sx q[0];
rz(-0.25704849) q[0];
rz(-2.2611179) q[1];
sx q[1];
rz(-0.6404666) q[1];
sx q[1];
rz(-1.2535198) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220164) q[0];
sx q[0];
rz(-1.5347693) q[0];
sx q[0];
rz(0.018112273) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80998151) q[2];
sx q[2];
rz(-1.4000386) q[2];
sx q[2];
rz(-2.8851938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5585564) q[1];
sx q[1];
rz(-0.7987403) q[1];
sx q[1];
rz(3.0786199) q[1];
x q[2];
rz(-1.569031) q[3];
sx q[3];
rz(-1.9179419) q[3];
sx q[3];
rz(-2.7070759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4132061) q[2];
sx q[2];
rz(-1.9325958) q[2];
sx q[2];
rz(-1.2665292) q[2];
rz(0.62147102) q[3];
sx q[3];
rz(-2.5340243) q[3];
sx q[3];
rz(2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40015873) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(-0.025064502) q[0];
rz(1.9301682) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(-2.8866344) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0023277442) q[0];
sx q[0];
rz(-0.049590913) q[0];
sx q[0];
rz(-3.0399486) q[0];
rz(2.0588875) q[2];
sx q[2];
rz(-2.1016663) q[2];
sx q[2];
rz(-1.7386029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.43422302) q[1];
sx q[1];
rz(-1.0138144) q[1];
sx q[1];
rz(-2.2051478) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9288152) q[3];
sx q[3];
rz(-2.1329974) q[3];
sx q[3];
rz(-0.24383814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24594626) q[2];
sx q[2];
rz(-1.0422372) q[2];
sx q[2];
rz(2.6825405) q[2];
rz(-1.6445271) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(-0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46886214) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(-3.0555308) q[0];
rz(-0.91880265) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(1.930621) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.092441) q[0];
sx q[0];
rz(-1.7730129) q[0];
sx q[0];
rz(-2.8114955) q[0];
x q[1];
rz(2.2077256) q[2];
sx q[2];
rz(-2.4709765) q[2];
sx q[2];
rz(-3.0702555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0748634) q[1];
sx q[1];
rz(-1.7742429) q[1];
sx q[1];
rz(-1.8977081) q[1];
rz(-pi) q[2];
rz(-1.9350697) q[3];
sx q[3];
rz(-1.2281872) q[3];
sx q[3];
rz(2.2628257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3226037) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(-2.4244579) q[2];
rz(2.1142193) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(-0.43924847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1502007) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(-0.014658654) q[0];
rz(-2.390059) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(1.470648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.178407) q[0];
sx q[0];
rz(-1.7589594) q[0];
sx q[0];
rz(1.7117731) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29701091) q[2];
sx q[2];
rz(-2.8737465) q[2];
sx q[2];
rz(-1.0849407) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5353706) q[1];
sx q[1];
rz(-0.41946966) q[1];
sx q[1];
rz(1.2910045) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5603546) q[3];
sx q[3];
rz(-0.34509788) q[3];
sx q[3];
rz(1.0853678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26059255) q[2];
sx q[2];
rz(-2.0286109) q[2];
sx q[2];
rz(-0.27349681) q[2];
rz(-1.1547487) q[3];
sx q[3];
rz(-0.65360779) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.038789373) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(-2.0661085) q[0];
rz(-3.0134046) q[1];
sx q[1];
rz(-0.85103858) q[1];
sx q[1];
rz(-2.7959965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1911666) q[0];
sx q[0];
rz(-1.7306402) q[0];
sx q[0];
rz(-2.2402083) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0021474) q[2];
sx q[2];
rz(-1.8446184) q[2];
sx q[2];
rz(1.3783) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3070655) q[1];
sx q[1];
rz(-2.5967) q[1];
sx q[1];
rz(0.12172555) q[1];
x q[2];
rz(-2.54353) q[3];
sx q[3];
rz(-1.4334213) q[3];
sx q[3];
rz(3.1118903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0547611) q[2];
sx q[2];
rz(-2.0589477) q[2];
sx q[2];
rz(1.6205988) q[2];
rz(-1.7704376) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(-2.7267406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3138251) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-0.23571043) q[0];
rz(2.56855) q[1];
sx q[1];
rz(-2.133281) q[1];
sx q[1];
rz(-1.7726353) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.154207) q[0];
sx q[0];
rz(-1.9999749) q[0];
sx q[0];
rz(-2.3230419) q[0];
rz(-0.70765453) q[2];
sx q[2];
rz(-0.81205149) q[2];
sx q[2];
rz(2.6132513) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6073709) q[1];
sx q[1];
rz(-2.1335619) q[1];
sx q[1];
rz(-2.9396179) q[1];
x q[2];
rz(-1.664613) q[3];
sx q[3];
rz(-1.929885) q[3];
sx q[3];
rz(2.1916568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3706751) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(0.053000432) q[2];
rz(-0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(2.2161765) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62405217) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(-1.7908295) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(2.9861957) q[2];
sx q[2];
rz(-0.70320319) q[2];
sx q[2];
rz(2.7413766) q[2];
rz(0.97196058) q[3];
sx q[3];
rz(-2.2836015) q[3];
sx q[3];
rz(-2.9290269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
