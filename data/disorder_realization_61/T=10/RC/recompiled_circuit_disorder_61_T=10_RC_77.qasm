OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4133889) q[0];
sx q[0];
rz(-1.1336741) q[0];
sx q[0];
rz(1.5925621) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(2.5974098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13574164) q[0];
sx q[0];
rz(-0.07677456) q[0];
sx q[0];
rz(2.6430921) q[0];
rz(-pi) q[1];
rz(2.7317023) q[2];
sx q[2];
rz(-3.0311243) q[2];
sx q[2];
rz(-2.8540749) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4197598) q[1];
sx q[1];
rz(-2.3082323) q[1];
sx q[1];
rz(-2.1268197) q[1];
x q[2];
rz(2.5892026) q[3];
sx q[3];
rz(-0.033621764) q[3];
sx q[3];
rz(-2.606483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0698174) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(1.3624181) q[2];
rz(-3.1133364) q[3];
sx q[3];
rz(-1.3794206) q[3];
sx q[3];
rz(-0.79022592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-0.70514482) q[0];
sx q[0];
rz(-1.171296) q[0];
rz(2.9303739) q[1];
sx q[1];
rz(-2.6995081) q[1];
sx q[1];
rz(1.404095) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728714) q[0];
sx q[0];
rz(-2.255548) q[0];
sx q[0];
rz(-1.0147592) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4046461) q[2];
sx q[2];
rz(-1.0937905) q[2];
sx q[2];
rz(1.1039066) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.055981936) q[1];
sx q[1];
rz(-1.8538875) q[1];
sx q[1];
rz(1.5468456) q[1];
rz(-1.7066129) q[3];
sx q[3];
rz(-1.0944301) q[3];
sx q[3];
rz(-1.7686896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30963787) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(-0.94397604) q[2];
rz(-0.55654636) q[3];
sx q[3];
rz(-2.6440933) q[3];
sx q[3];
rz(-1.336162) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29207644) q[0];
sx q[0];
rz(-0.76382604) q[0];
sx q[0];
rz(1.4235494) q[0];
rz(-2.3220093) q[1];
sx q[1];
rz(-0.81033605) q[1];
sx q[1];
rz(0.56366411) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1467928) q[0];
sx q[0];
rz(-2.1437862) q[0];
sx q[0];
rz(1.6550199) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7965505) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(3.0622481) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1395531) q[1];
sx q[1];
rz(-2.4783652) q[1];
sx q[1];
rz(1.0242277) q[1];
rz(1.1590957) q[3];
sx q[3];
rz(-0.99038306) q[3];
sx q[3];
rz(-1.8713649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6391969) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(2.0813265) q[2];
rz(0.075573102) q[3];
sx q[3];
rz(-0.85791701) q[3];
sx q[3];
rz(-3.0269572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.32325) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(0.19317214) q[0];
rz(1.5441783) q[1];
sx q[1];
rz(-1.1491821) q[1];
sx q[1];
rz(0.65778041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5415216) q[0];
sx q[0];
rz(-2.9642448) q[0];
sx q[0];
rz(2.970201) q[0];
rz(-pi) q[1];
rz(-2.0006318) q[2];
sx q[2];
rz(-1.8550711) q[2];
sx q[2];
rz(-0.62361275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.12280497) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(2.8469574) q[1];
rz(-pi) q[2];
rz(1.863402) q[3];
sx q[3];
rz(-1.3232104) q[3];
sx q[3];
rz(2.2729371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0458935) q[2];
sx q[2];
rz(-2.7313488) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(2.0783157) q[3];
sx q[3];
rz(-1.9789109) q[3];
sx q[3];
rz(-1.8803546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3445774) q[0];
sx q[0];
rz(-0.319096) q[0];
sx q[0];
rz(-2.1176594) q[0];
rz(0.78760415) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(0.62686282) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.244827) q[0];
sx q[0];
rz(-1.5540431) q[0];
sx q[0];
rz(1.6003952) q[0];
rz(-0.31693964) q[2];
sx q[2];
rz(-1.8570515) q[2];
sx q[2];
rz(1.1405108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8265241) q[1];
sx q[1];
rz(-1.2722004) q[1];
sx q[1];
rz(-1.8930757) q[1];
rz(-pi) q[2];
rz(0.12368006) q[3];
sx q[3];
rz(-1.1808174) q[3];
sx q[3];
rz(2.9434413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2287067) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(2.1234925) q[2];
rz(-2.5300238) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6927032) q[0];
sx q[0];
rz(-1.3494116) q[0];
sx q[0];
rz(-1.942379) q[0];
rz(-1.1692283) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.823846) q[0];
sx q[0];
rz(-2.4577603) q[0];
sx q[0];
rz(-0.86156396) q[0];
rz(-pi) q[1];
rz(2.3338823) q[2];
sx q[2];
rz(-2.0732023) q[2];
sx q[2];
rz(-3.0903357) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2734087) q[1];
sx q[1];
rz(-1.5460099) q[1];
sx q[1];
rz(-2.1048057) q[1];
rz(-1.2732029) q[3];
sx q[3];
rz(-1.9801567) q[3];
sx q[3];
rz(0.62664947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4679608) q[2];
sx q[2];
rz(-1.3053852) q[2];
sx q[2];
rz(-1.2754053) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4642898) q[0];
sx q[0];
rz(-1.3339366) q[0];
sx q[0];
rz(-1.4087079) q[0];
rz(0.45577058) q[1];
sx q[1];
rz(-0.20142889) q[1];
sx q[1];
rz(-1.2021525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81284886) q[0];
sx q[0];
rz(-2.3658731) q[0];
sx q[0];
rz(-2.1562955) q[0];
rz(-0.41263327) q[2];
sx q[2];
rz(-2.3473783) q[2];
sx q[2];
rz(-2.7318294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.251293) q[1];
sx q[1];
rz(-1.4874465) q[1];
sx q[1];
rz(3.0492196) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2644516) q[3];
sx q[3];
rz(-1.3925941) q[3];
sx q[3];
rz(-1.7314535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.000164) q[2];
sx q[2];
rz(-1.8464073) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(-2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(-1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0705567) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-1.113168) q[0];
rz(0.69560266) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.5323458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5672011) q[0];
sx q[0];
rz(-2.6041457) q[0];
sx q[0];
rz(0.48899942) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74560994) q[2];
sx q[2];
rz(-0.62586212) q[2];
sx q[2];
rz(-1.2127884) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.53837) q[1];
sx q[1];
rz(-2.6909628) q[1];
sx q[1];
rz(1.6166812) q[1];
rz(1.6892151) q[3];
sx q[3];
rz(-0.97202557) q[3];
sx q[3];
rz(1.0373725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1476851) q[2];
sx q[2];
rz(-0.2299749) q[2];
sx q[2];
rz(-0.85419401) q[2];
rz(-1.2285852) q[3];
sx q[3];
rz(-1.5766141) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937186) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(-1.9539072) q[1];
sx q[1];
rz(-0.70078754) q[1];
sx q[1];
rz(1.2845576) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4391543) q[0];
sx q[0];
rz(-2.2674773) q[0];
sx q[0];
rz(-0.37581635) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8753488) q[2];
sx q[2];
rz(-1.9387445) q[2];
sx q[2];
rz(3.0014696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6397275) q[1];
sx q[1];
rz(-1.1266536) q[1];
sx q[1];
rz(0.72413866) q[1];
rz(-pi) q[2];
rz(-1.5981538) q[3];
sx q[3];
rz(-1.9606855) q[3];
sx q[3];
rz(-1.1297806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5294042) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(2.5908296) q[2];
rz(-0.70358706) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(-1.6350869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4187014) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(0.044145949) q[0];
rz(1.6126397) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(0.77967656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44256193) q[0];
sx q[0];
rz(-0.5578707) q[0];
sx q[0];
rz(-2.872422) q[0];
x q[1];
rz(-2.6780307) q[2];
sx q[2];
rz(-2.2710685) q[2];
sx q[2];
rz(-2.314032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.838678) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.7811437) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3603504) q[3];
sx q[3];
rz(-1.9133854) q[3];
sx q[3];
rz(-2.4398746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49446517) q[2];
sx q[2];
rz(-2.4202042) q[2];
sx q[2];
rz(1.5987827) q[2];
rz(2.4370082) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(0.012044756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71173944) q[0];
sx q[0];
rz(-1.585351) q[0];
sx q[0];
rz(2.4899695) q[0];
rz(-1.343887) q[1];
sx q[1];
rz(-1.4812891) q[1];
sx q[1];
rz(-0.67566009) q[1];
rz(-0.96991878) q[2];
sx q[2];
rz(-2.6261332) q[2];
sx q[2];
rz(0.020608227) q[2];
rz(-1.4383437) q[3];
sx q[3];
rz(-2.2046897) q[3];
sx q[3];
rz(-0.71837367) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];