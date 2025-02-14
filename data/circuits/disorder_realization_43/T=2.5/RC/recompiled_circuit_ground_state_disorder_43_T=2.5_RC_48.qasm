OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3620152) q[0];
sx q[0];
rz(-2.9147122) q[0];
sx q[0];
rz(1.3751295) q[0];
rz(3.0472164) q[1];
sx q[1];
rz(-0.91369319) q[1];
sx q[1];
rz(0.28092608) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0931647) q[0];
sx q[0];
rz(-2.4602088) q[0];
sx q[0];
rz(1.3171893) q[0];
rz(-pi) q[1];
rz(-1.2035349) q[2];
sx q[2];
rz(-0.31031552) q[2];
sx q[2];
rz(0.65904891) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0795796) q[1];
sx q[1];
rz(-1.8342736) q[1];
sx q[1];
rz(2.0017712) q[1];
x q[2];
rz(-0.51961629) q[3];
sx q[3];
rz(-0.63881385) q[3];
sx q[3];
rz(0.45376247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1137587) q[2];
sx q[2];
rz(-1.4049302) q[2];
sx q[2];
rz(-3.0244381) q[2];
rz(-2.8095918) q[3];
sx q[3];
rz(-0.74688512) q[3];
sx q[3];
rz(2.3565256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4431045) q[0];
sx q[0];
rz(-2.3790058) q[0];
sx q[0];
rz(2.7681328) q[0];
rz(-2.7442878) q[1];
sx q[1];
rz(-1.1445069) q[1];
sx q[1];
rz(-0.73748803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0076548) q[0];
sx q[0];
rz(-2.2077401) q[0];
sx q[0];
rz(-2.4600814) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5869484) q[2];
sx q[2];
rz(-2.0337542) q[2];
sx q[2];
rz(2.2210768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7389679) q[1];
sx q[1];
rz(-2.1161656) q[1];
sx q[1];
rz(1.7508372) q[1];
x q[2];
rz(-1.9225538) q[3];
sx q[3];
rz(-2.2969807) q[3];
sx q[3];
rz(1.6399469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.069933683) q[2];
sx q[2];
rz(-1.5295014) q[2];
sx q[2];
rz(2.7640479) q[2];
rz(-0.30970445) q[3];
sx q[3];
rz(-2.0524502) q[3];
sx q[3];
rz(-1.6449876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6250703) q[0];
sx q[0];
rz(-0.83928883) q[0];
sx q[0];
rz(1.2155493) q[0];
rz(-2.1274321) q[1];
sx q[1];
rz(-1.4233669) q[1];
sx q[1];
rz(-0.97022143) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6785203) q[0];
sx q[0];
rz(-2.6404535) q[0];
sx q[0];
rz(1.4864413) q[0];
rz(0.54527905) q[2];
sx q[2];
rz(-2.1542794) q[2];
sx q[2];
rz(0.80697053) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41170317) q[1];
sx q[1];
rz(-1.7474993) q[1];
sx q[1];
rz(-1.9167625) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30556314) q[3];
sx q[3];
rz(-0.36763469) q[3];
sx q[3];
rz(-0.43722269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0936475) q[2];
sx q[2];
rz(-2.2990172) q[2];
sx q[2];
rz(-2.688431) q[2];
rz(-1.9122745) q[3];
sx q[3];
rz(-1.2776351) q[3];
sx q[3];
rz(2.9696828) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1969084) q[0];
sx q[0];
rz(-0.6854282) q[0];
sx q[0];
rz(0.36566439) q[0];
rz(0.91224313) q[1];
sx q[1];
rz(-1.279) q[1];
sx q[1];
rz(-2.7361187) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.643754) q[0];
sx q[0];
rz(-0.23786834) q[0];
sx q[0];
rz(-2.0980623) q[0];
rz(-pi) q[1];
rz(2.916476) q[2];
sx q[2];
rz(-0.76864132) q[2];
sx q[2];
rz(-0.74736881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7995389) q[1];
sx q[1];
rz(-2.5785682) q[1];
sx q[1];
rz(2.8232226) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41678352) q[3];
sx q[3];
rz(-2.458771) q[3];
sx q[3];
rz(0.58931755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0696062) q[2];
sx q[2];
rz(-1.7013841) q[2];
sx q[2];
rz(-1.5427422) q[2];
rz(0.37374464) q[3];
sx q[3];
rz(-1.475324) q[3];
sx q[3];
rz(-1.4603978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54295802) q[0];
sx q[0];
rz(-0.77474189) q[0];
sx q[0];
rz(-2.4454818) q[0];
rz(-1.8822582) q[1];
sx q[1];
rz(-2.0399317) q[1];
sx q[1];
rz(-2.7764244) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0232622) q[0];
sx q[0];
rz(-2.0666762) q[0];
sx q[0];
rz(2.0642679) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1579726) q[2];
sx q[2];
rz(-1.8637805) q[2];
sx q[2];
rz(-2.3940046) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3427178) q[1];
sx q[1];
rz(-1.5694071) q[1];
sx q[1];
rz(0.54552127) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5138084) q[3];
sx q[3];
rz(-1.2300228) q[3];
sx q[3];
rz(-0.22751156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.460707) q[2];
sx q[2];
rz(-1.1628217) q[2];
sx q[2];
rz(-2.9648901) q[2];
rz(1.3867311) q[3];
sx q[3];
rz(-2.0818129) q[3];
sx q[3];
rz(-2.7669014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52727592) q[0];
sx q[0];
rz(-2.2820331) q[0];
sx q[0];
rz(1.5981307) q[0];
rz(1.3044926) q[1];
sx q[1];
rz(-1.0544798) q[1];
sx q[1];
rz(-0.80610448) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4325754) q[0];
sx q[0];
rz(-0.78993778) q[0];
sx q[0];
rz(-2.8054601) q[0];
rz(0.51050425) q[2];
sx q[2];
rz(-2.1091828) q[2];
sx q[2];
rz(-1.5988072) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.66530217) q[1];
sx q[1];
rz(-0.21033439) q[1];
sx q[1];
rz(1.8108941) q[1];
rz(-pi) q[2];
rz(-1.0433572) q[3];
sx q[3];
rz(-1.0602078) q[3];
sx q[3];
rz(3.1289738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.97258893) q[2];
sx q[2];
rz(-2.6469595) q[2];
sx q[2];
rz(0.063610323) q[2];
rz(-0.58498597) q[3];
sx q[3];
rz(-1.4002742) q[3];
sx q[3];
rz(-1.5144279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1503898) q[0];
sx q[0];
rz(-1.7504033) q[0];
sx q[0];
rz(-0.72917953) q[0];
rz(-1.179262) q[1];
sx q[1];
rz(-2.3919892) q[1];
sx q[1];
rz(2.8470305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33074327) q[0];
sx q[0];
rz(-2.0669575) q[0];
sx q[0];
rz(-1.2788354) q[0];
rz(-pi) q[1];
rz(1.137758) q[2];
sx q[2];
rz(-1.6707509) q[2];
sx q[2];
rz(-1.86509) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.197261) q[1];
sx q[1];
rz(-0.75350584) q[1];
sx q[1];
rz(1.7427518) q[1];
rz(1.1001281) q[3];
sx q[3];
rz(-2.4146955) q[3];
sx q[3];
rz(-1.5125546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4906759) q[2];
sx q[2];
rz(-1.4399521) q[2];
sx q[2];
rz(0.73224625) q[2];
rz(-2.2722774) q[3];
sx q[3];
rz(-1.1704051) q[3];
sx q[3];
rz(1.4066345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9100087) q[0];
sx q[0];
rz(-1.5861479) q[0];
sx q[0];
rz(0.79767942) q[0];
rz(-0.32294598) q[1];
sx q[1];
rz(-0.99635092) q[1];
sx q[1];
rz(1.4083883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7777611) q[0];
sx q[0];
rz(-2.1182502) q[0];
sx q[0];
rz(1.8108032) q[0];
rz(-2.0261835) q[2];
sx q[2];
rz(-0.50297996) q[2];
sx q[2];
rz(3.0969381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6646386) q[1];
sx q[1];
rz(-0.79915291) q[1];
sx q[1];
rz(0.33070143) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3207758) q[3];
sx q[3];
rz(-1.3969706) q[3];
sx q[3];
rz(2.6717693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8352167) q[2];
sx q[2];
rz(-1.110346) q[2];
sx q[2];
rz(0.19732538) q[2];
rz(-0.80162445) q[3];
sx q[3];
rz(-2.1620731) q[3];
sx q[3];
rz(-0.42158034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9824958) q[0];
sx q[0];
rz(-1.4747341) q[0];
sx q[0];
rz(0.91484797) q[0];
rz(-0.21028701) q[1];
sx q[1];
rz(-0.57840127) q[1];
sx q[1];
rz(-2.4519144) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1062968) q[0];
sx q[0];
rz(-1.4778504) q[0];
sx q[0];
rz(-0.22648099) q[0];
rz(0.014691512) q[2];
sx q[2];
rz(-1.6882875) q[2];
sx q[2];
rz(-1.9368287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2336809) q[1];
sx q[1];
rz(-2.4123976) q[1];
sx q[1];
rz(-2.3883874) q[1];
rz(-pi) q[2];
rz(2.0965936) q[3];
sx q[3];
rz(-1.8711578) q[3];
sx q[3];
rz(-1.5142358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0939533) q[2];
sx q[2];
rz(-0.84838212) q[2];
sx q[2];
rz(0.32996714) q[2];
rz(1.0432358) q[3];
sx q[3];
rz(-1.7236575) q[3];
sx q[3];
rz(-0.51658336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.030815) q[0];
sx q[0];
rz(-1.8106221) q[0];
sx q[0];
rz(-1.0571085) q[0];
rz(-0.2541751) q[1];
sx q[1];
rz(-1.1421685) q[1];
sx q[1];
rz(-2.4370297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6374311) q[0];
sx q[0];
rz(-1.1011718) q[0];
sx q[0];
rz(-0.35029098) q[0];
rz(-1.754934) q[2];
sx q[2];
rz(-1.6312851) q[2];
sx q[2];
rz(0.16512251) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8125638) q[1];
sx q[1];
rz(-2.6575251) q[1];
sx q[1];
rz(-1.1776393) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4289565) q[3];
sx q[3];
rz(-2.8691926) q[3];
sx q[3];
rz(-0.89100641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5083984) q[2];
sx q[2];
rz(-1.1722379) q[2];
sx q[2];
rz(-2.634826) q[2];
rz(0.1861598) q[3];
sx q[3];
rz(-0.23701826) q[3];
sx q[3];
rz(2.6059634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7645466) q[0];
sx q[0];
rz(-1.4519539) q[0];
sx q[0];
rz(-2.0013381) q[0];
rz(-2.2346732) q[1];
sx q[1];
rz(-2.4005371) q[1];
sx q[1];
rz(-1.3225318) q[1];
rz(2.2494153) q[2];
sx q[2];
rz(-1.1363251) q[2];
sx q[2];
rz(-2.688495) q[2];
rz(-1.3033397) q[3];
sx q[3];
rz(-1.8831913) q[3];
sx q[3];
rz(2.5358653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
