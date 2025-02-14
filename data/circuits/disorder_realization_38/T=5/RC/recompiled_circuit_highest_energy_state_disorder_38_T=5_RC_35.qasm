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
rz(0.44904798) q[0];
sx q[0];
rz(4.2807978) q[0];
sx q[0];
rz(9.8674404) q[0];
rz(-2.9709385) q[1];
sx q[1];
rz(-0.79167384) q[1];
sx q[1];
rz(-0.72905529) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9937482) q[0];
sx q[0];
rz(-1.1708852) q[0];
sx q[0];
rz(-0.73109674) q[0];
rz(-pi) q[1];
rz(-1.5716971) q[2];
sx q[2];
rz(-1.9605165) q[2];
sx q[2];
rz(-2.2665521) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3786582) q[1];
sx q[1];
rz(-2.4097917) q[1];
sx q[1];
rz(-1.9373489) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40790848) q[3];
sx q[3];
rz(-2.1589557) q[3];
sx q[3];
rz(-1.5666636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.86058229) q[2];
sx q[2];
rz(-3.0594825) q[2];
sx q[2];
rz(2.3561467) q[2];
rz(-0.9663409) q[3];
sx q[3];
rz(-2.0735425) q[3];
sx q[3];
rz(2.5016224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53168374) q[0];
sx q[0];
rz(-0.56782472) q[0];
sx q[0];
rz(0.63572836) q[0];
rz(2.5081778) q[1];
sx q[1];
rz(-0.76786357) q[1];
sx q[1];
rz(0.33087081) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61900789) q[0];
sx q[0];
rz(-0.96946335) q[0];
sx q[0];
rz(1.1669789) q[0];
x q[1];
rz(1.6236247) q[2];
sx q[2];
rz(-0.46762782) q[2];
sx q[2];
rz(-2.8421081) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91869044) q[1];
sx q[1];
rz(-1.6044334) q[1];
sx q[1];
rz(-1.6292389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5583573) q[3];
sx q[3];
rz(-1.6553989) q[3];
sx q[3];
rz(-2.0432664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6364381) q[2];
sx q[2];
rz(-2.324489) q[2];
sx q[2];
rz(2.6804697) q[2];
rz(0.72502208) q[3];
sx q[3];
rz(-2.6889763) q[3];
sx q[3];
rz(-2.4729474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4918936) q[0];
sx q[0];
rz(-1.7576341) q[0];
sx q[0];
rz(2.6032676) q[0];
rz(3.1321943) q[1];
sx q[1];
rz(-0.63176578) q[1];
sx q[1];
rz(1.1993154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56012404) q[0];
sx q[0];
rz(-2.631979) q[0];
sx q[0];
rz(1.4790003) q[0];
rz(2.968987) q[2];
sx q[2];
rz(-0.69059082) q[2];
sx q[2];
rz(-2.0334854) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4863094) q[1];
sx q[1];
rz(-1.9260672) q[1];
sx q[1];
rz(2.1706697) q[1];
rz(-pi) q[2];
rz(1.2869419) q[3];
sx q[3];
rz(-2.3945269) q[3];
sx q[3];
rz(3.0892854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37375307) q[2];
sx q[2];
rz(-1.491863) q[2];
sx q[2];
rz(-0.72244942) q[2];
rz(0.59822285) q[3];
sx q[3];
rz(-3.1066419) q[3];
sx q[3];
rz(1.2178864) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55906051) q[0];
sx q[0];
rz(-1.8836972) q[0];
sx q[0];
rz(3.1225358) q[0];
rz(1.4590774) q[1];
sx q[1];
rz(-0.2225114) q[1];
sx q[1];
rz(2.6313475) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6551901) q[0];
sx q[0];
rz(-2.1684747) q[0];
sx q[0];
rz(-2.8199353) q[0];
rz(2.4195477) q[2];
sx q[2];
rz(-1.8068442) q[2];
sx q[2];
rz(-1.2942734) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7170808) q[1];
sx q[1];
rz(-1.9766909) q[1];
sx q[1];
rz(0.67274977) q[1];
x q[2];
rz(1.7137547) q[3];
sx q[3];
rz(-1.1605439) q[3];
sx q[3];
rz(-2.9825135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14429188) q[2];
sx q[2];
rz(-1.68953) q[2];
sx q[2];
rz(-0.22225456) q[2];
rz(3.0620767) q[3];
sx q[3];
rz(-0.54602081) q[3];
sx q[3];
rz(-2.5047746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6154489) q[0];
sx q[0];
rz(-2.0455102) q[0];
sx q[0];
rz(1.8002321) q[0];
rz(0.51097956) q[1];
sx q[1];
rz(-1.363441) q[1];
sx q[1];
rz(2.708639) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6658185) q[0];
sx q[0];
rz(-1.6278212) q[0];
sx q[0];
rz(0.25754575) q[0];
rz(1.9741945) q[2];
sx q[2];
rz(-1.7678542) q[2];
sx q[2];
rz(2.6464406) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4652327) q[1];
sx q[1];
rz(-1.3238878) q[1];
sx q[1];
rz(1.2311155) q[1];
x q[2];
rz(-0.28750395) q[3];
sx q[3];
rz(-2.3477481) q[3];
sx q[3];
rz(-0.35121976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13010919) q[2];
sx q[2];
rz(-2.2138962) q[2];
sx q[2];
rz(-0.9978655) q[2];
rz(0.70895553) q[3];
sx q[3];
rz(-0.38781375) q[3];
sx q[3];
rz(-0.97878218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82584941) q[0];
sx q[0];
rz(-2.5542927) q[0];
sx q[0];
rz(-0.89383268) q[0];
rz(-1.029344) q[1];
sx q[1];
rz(-2.4004332) q[1];
sx q[1];
rz(3.0267874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66291565) q[0];
sx q[0];
rz(-1.3548242) q[0];
sx q[0];
rz(1.8195282) q[0];
rz(-pi) q[1];
rz(-0.46415879) q[2];
sx q[2];
rz(-0.69550976) q[2];
sx q[2];
rz(1.1087195) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8821174) q[1];
sx q[1];
rz(-2.7434556) q[1];
sx q[1];
rz(-0.78146259) q[1];
rz(1.6345665) q[3];
sx q[3];
rz(-1.290451) q[3];
sx q[3];
rz(0.75322039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.43742314) q[2];
sx q[2];
rz(-2.2544474) q[2];
sx q[2];
rz(-2.9278921) q[2];
rz(-2.5050488) q[3];
sx q[3];
rz(-0.53037363) q[3];
sx q[3];
rz(-0.73591939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(0.32695025) q[0];
sx q[0];
rz(-2.3982168) q[0];
sx q[0];
rz(3.0315234) q[0];
rz(-3.0649109) q[1];
sx q[1];
rz(-2.2814467) q[1];
sx q[1];
rz(-1.1383879) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6784043) q[0];
sx q[0];
rz(-2.7240755) q[0];
sx q[0];
rz(-0.56708401) q[0];
rz(0.78779548) q[2];
sx q[2];
rz(-1.5681364) q[2];
sx q[2];
rz(-0.12622368) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7775007) q[1];
sx q[1];
rz(-2.052825) q[1];
sx q[1];
rz(-2.9338323) q[1];
x q[2];
rz(-2.4289598) q[3];
sx q[3];
rz(-0.32100916) q[3];
sx q[3];
rz(2.6965303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8818605) q[2];
sx q[2];
rz(-2.0596108) q[2];
sx q[2];
rz(-2.760375) q[2];
rz(-3.0242053) q[3];
sx q[3];
rz(-2.6365247) q[3];
sx q[3];
rz(-0.086732619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7748902) q[0];
sx q[0];
rz(-2.369304) q[0];
sx q[0];
rz(0.49852398) q[0];
rz(0.82930928) q[1];
sx q[1];
rz(-0.81559759) q[1];
sx q[1];
rz(-0.097749762) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691188) q[0];
sx q[0];
rz(-2.0017255) q[0];
sx q[0];
rz(2.8043967) q[0];
rz(-1.421861) q[2];
sx q[2];
rz(-1.2342936) q[2];
sx q[2];
rz(0.83068427) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4227998) q[1];
sx q[1];
rz(-0.96578078) q[1];
sx q[1];
rz(2.2284719) q[1];
x q[2];
rz(0.0093383647) q[3];
sx q[3];
rz(-1.2492325) q[3];
sx q[3];
rz(-0.044794785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30663124) q[2];
sx q[2];
rz(-1.9913048) q[2];
sx q[2];
rz(-0.59070307) q[2];
rz(0.090204209) q[3];
sx q[3];
rz(-0.43961757) q[3];
sx q[3];
rz(-2.2265767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0331994) q[0];
sx q[0];
rz(-2.2985701) q[0];
sx q[0];
rz(0.75476187) q[0];
rz(-2.8056878) q[1];
sx q[1];
rz(-2.5867808) q[1];
sx q[1];
rz(1.0364484) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51506804) q[0];
sx q[0];
rz(-1.4534411) q[0];
sx q[0];
rz(2.4548454) q[0];
rz(-pi) q[1];
rz(-0.48123863) q[2];
sx q[2];
rz(-1.5889152) q[2];
sx q[2];
rz(2.6157524) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8120665) q[1];
sx q[1];
rz(-0.5181075) q[1];
sx q[1];
rz(0.078703367) q[1];
rz(-pi) q[2];
x q[2];
rz(0.01364661) q[3];
sx q[3];
rz(-1.3411152) q[3];
sx q[3];
rz(-1.0245293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0466517) q[2];
sx q[2];
rz(-1.4163821) q[2];
sx q[2];
rz(2.6851192) q[2];
rz(-2.5392695) q[3];
sx q[3];
rz(-2.9538302) q[3];
sx q[3];
rz(-2.203234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8896821) q[0];
sx q[0];
rz(-1.849181) q[0];
sx q[0];
rz(-0.46345261) q[0];
rz(-3.1262596) q[1];
sx q[1];
rz(-0.066611193) q[1];
sx q[1];
rz(2.8212246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.244359) q[0];
sx q[0];
rz(-0.16288745) q[0];
sx q[0];
rz(-2.6475859) q[0];
rz(-1.8202844) q[2];
sx q[2];
rz(-1.87566) q[2];
sx q[2];
rz(2.9314205) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3829455) q[1];
sx q[1];
rz(-1.0898814) q[1];
sx q[1];
rz(0.8713027) q[1];
x q[2];
rz(1.2804561) q[3];
sx q[3];
rz(-1.2429894) q[3];
sx q[3];
rz(2.7442055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1128803) q[2];
sx q[2];
rz(-0.72673231) q[2];
sx q[2];
rz(1.0090562) q[2];
rz(0.79695898) q[3];
sx q[3];
rz(-1.8877441) q[3];
sx q[3];
rz(-2.8033281) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0924031) q[0];
sx q[0];
rz(-1.5143464) q[0];
sx q[0];
rz(-1.2595246) q[0];
rz(-3.0567723) q[1];
sx q[1];
rz(-1.6592574) q[1];
sx q[1];
rz(1.4316373) q[1];
rz(1.8126429) q[2];
sx q[2];
rz(-1.6326661) q[2];
sx q[2];
rz(2.2018955) q[2];
rz(-2.8390213) q[3];
sx q[3];
rz(-1.745001) q[3];
sx q[3];
rz(0.59636084) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
