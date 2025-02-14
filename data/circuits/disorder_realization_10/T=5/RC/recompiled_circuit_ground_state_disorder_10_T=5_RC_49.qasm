OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4781696) q[0];
sx q[0];
rz(-0.91376758) q[0];
sx q[0];
rz(-2.0765641) q[0];
rz(1.6051259) q[1];
sx q[1];
rz(-1.1054339) q[1];
sx q[1];
rz(0.096535834) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6632412) q[0];
sx q[0];
rz(-2.6697201) q[0];
sx q[0];
rz(-2.0357516) q[0];
rz(-0.41806721) q[2];
sx q[2];
rz(-2.5804248) q[2];
sx q[2];
rz(0.18218606) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9241544) q[1];
sx q[1];
rz(-0.70637843) q[1];
sx q[1];
rz(-0.46831727) q[1];
rz(-1.4289342) q[3];
sx q[3];
rz(-0.58462287) q[3];
sx q[3];
rz(-1.4277924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1221293) q[2];
sx q[2];
rz(-2.9563642) q[2];
sx q[2];
rz(1.2629925) q[2];
rz(-0.39877912) q[3];
sx q[3];
rz(-1.140927) q[3];
sx q[3];
rz(-2.1122475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(1.0313834) q[0];
sx q[0];
rz(-0.6321913) q[0];
sx q[0];
rz(-1.6803918) q[0];
rz(-3.0714495) q[1];
sx q[1];
rz(-1.5927916) q[1];
sx q[1];
rz(-0.4313012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3628497) q[0];
sx q[0];
rz(-1.5093137) q[0];
sx q[0];
rz(-3.132122) q[0];
rz(-pi) q[1];
rz(-1.9877983) q[2];
sx q[2];
rz(-1.2975311) q[2];
sx q[2];
rz(2.5962256) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97535061) q[1];
sx q[1];
rz(-2.0845452) q[1];
sx q[1];
rz(2.3523037) q[1];
rz(-0.43092315) q[3];
sx q[3];
rz(-0.049073372) q[3];
sx q[3];
rz(1.2653093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.18616072) q[2];
sx q[2];
rz(-1.4107076) q[2];
sx q[2];
rz(0.53039941) q[2];
rz(2.0576599) q[3];
sx q[3];
rz(-1.089774) q[3];
sx q[3];
rz(-1.2497905) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1029866) q[0];
sx q[0];
rz(-2.5654721) q[0];
sx q[0];
rz(1.2023793) q[0];
rz(-1.0105969) q[1];
sx q[1];
rz(-1.8162411) q[1];
sx q[1];
rz(0.030166322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9011616) q[0];
sx q[0];
rz(-1.636292) q[0];
sx q[0];
rz(-1.5237048) q[0];
x q[1];
rz(-3.1254076) q[2];
sx q[2];
rz(-1.5183518) q[2];
sx q[2];
rz(-1.7982782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2304014) q[1];
sx q[1];
rz(-1.4433858) q[1];
sx q[1];
rz(-0.34613737) q[1];
rz(-pi) q[2];
x q[2];
rz(2.143489) q[3];
sx q[3];
rz(-1.4855322) q[3];
sx q[3];
rz(-2.5424438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5304337) q[2];
sx q[2];
rz(-2.9035089) q[2];
sx q[2];
rz(-0.37453026) q[2];
rz(-2.020906) q[3];
sx q[3];
rz(-1.4846385) q[3];
sx q[3];
rz(-0.69168958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.94367868) q[0];
sx q[0];
rz(-1.2926084) q[0];
sx q[0];
rz(1.9669272) q[0];
rz(1.8265751) q[1];
sx q[1];
rz(-2.0349793) q[1];
sx q[1];
rz(-0.16474251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667906) q[0];
sx q[0];
rz(-2.1336477) q[0];
sx q[0];
rz(0.99250162) q[0];
rz(1.7354749) q[2];
sx q[2];
rz(-1.0633755) q[2];
sx q[2];
rz(-0.31854445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81329423) q[1];
sx q[1];
rz(-0.61200145) q[1];
sx q[1];
rz(-2.185142) q[1];
x q[2];
rz(1.4993787) q[3];
sx q[3];
rz(-1.6096228) q[3];
sx q[3];
rz(-3.0049472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.69079196) q[2];
sx q[2];
rz(-1.5361293) q[2];
sx q[2];
rz(-0.50626051) q[2];
rz(-2.6311503) q[3];
sx q[3];
rz(-2.3531395) q[3];
sx q[3];
rz(-2.3300664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4617758) q[0];
sx q[0];
rz(-2.8051069) q[0];
sx q[0];
rz(0.80999723) q[0];
rz(3.0432155) q[1];
sx q[1];
rz(-2.7521303) q[1];
sx q[1];
rz(-0.32442763) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3220935) q[0];
sx q[0];
rz(-3.0190912) q[0];
sx q[0];
rz(1.5550645) q[0];
x q[1];
rz(-2.5971707) q[2];
sx q[2];
rz(-2.8834497) q[2];
sx q[2];
rz(-2.7552257) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4091332) q[1];
sx q[1];
rz(-0.9995102) q[1];
sx q[1];
rz(1.3582699) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3611322) q[3];
sx q[3];
rz(-2.0381322) q[3];
sx q[3];
rz(-1.0059822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.35679945) q[2];
sx q[2];
rz(-1.5609799) q[2];
sx q[2];
rz(0.55850935) q[2];
rz(-1.8465346) q[3];
sx q[3];
rz(-1.4174771) q[3];
sx q[3];
rz(0.3524802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2297939) q[0];
sx q[0];
rz(-0.42775446) q[0];
sx q[0];
rz(-1.7708923) q[0];
rz(1.7478583) q[1];
sx q[1];
rz(-1.0153898) q[1];
sx q[1];
rz(0.14399354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27731652) q[0];
sx q[0];
rz(-1.6747518) q[0];
sx q[0];
rz(-0.070965537) q[0];
x q[1];
rz(2.0259816) q[2];
sx q[2];
rz(-1.7894288) q[2];
sx q[2];
rz(-1.3022547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.86477937) q[1];
sx q[1];
rz(-2.1739156) q[1];
sx q[1];
rz(3.0052207) q[1];
rz(-pi) q[2];
rz(2.6089588) q[3];
sx q[3];
rz(-0.71511474) q[3];
sx q[3];
rz(-1.2330513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7250942) q[2];
sx q[2];
rz(-2.7855253) q[2];
sx q[2];
rz(-2.1605675) q[2];
rz(-0.52870005) q[3];
sx q[3];
rz(-1.3542391) q[3];
sx q[3];
rz(2.5813812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8037146) q[0];
sx q[0];
rz(-2.3382472) q[0];
sx q[0];
rz(-1.9377608) q[0];
rz(-0.0040668049) q[1];
sx q[1];
rz(-1.4264359) q[1];
sx q[1];
rz(0.34122658) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.041164) q[0];
sx q[0];
rz(-0.89027864) q[0];
sx q[0];
rz(-2.1245703) q[0];
x q[1];
rz(-2.262142) q[2];
sx q[2];
rz(-1.571889) q[2];
sx q[2];
rz(0.06070965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3918889) q[1];
sx q[1];
rz(-0.14924696) q[1];
sx q[1];
rz(0.061953739) q[1];
rz(2.7988966) q[3];
sx q[3];
rz(-1.7083941) q[3];
sx q[3];
rz(-2.9901341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61871201) q[2];
sx q[2];
rz(-1.2423542) q[2];
sx q[2];
rz(-0.27113554) q[2];
rz(1.5118269) q[3];
sx q[3];
rz(-0.96704331) q[3];
sx q[3];
rz(0.40201521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4790344) q[0];
sx q[0];
rz(-0.56773487) q[0];
sx q[0];
rz(0.16831368) q[0];
rz(1.0100826) q[1];
sx q[1];
rz(-1.964566) q[1];
sx q[1];
rz(-3.041306) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0675285) q[0];
sx q[0];
rz(-2.7180928) q[0];
sx q[0];
rz(0.74396043) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9434154) q[2];
sx q[2];
rz(-0.74286425) q[2];
sx q[2];
rz(-2.4695549) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.058096407) q[1];
sx q[1];
rz(-1.5952949) q[1];
sx q[1];
rz(-0.20170881) q[1];
rz(-2.4504205) q[3];
sx q[3];
rz(-0.6414957) q[3];
sx q[3];
rz(-0.79090276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4526796) q[2];
sx q[2];
rz(-2.6426688) q[2];
sx q[2];
rz(-1.4208687) q[2];
rz(1.0639327) q[3];
sx q[3];
rz(-1.7199793) q[3];
sx q[3];
rz(3.0498144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7559779) q[0];
sx q[0];
rz(-1.6843963) q[0];
sx q[0];
rz(-3.1047367) q[0];
rz(-1.256975) q[1];
sx q[1];
rz(-1.5024065) q[1];
sx q[1];
rz(-1.3704376) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12808558) q[0];
sx q[0];
rz(-2.8261746) q[0];
sx q[0];
rz(-1.0491788) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37439926) q[2];
sx q[2];
rz(-2.7485195) q[2];
sx q[2];
rz(-3.0469325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3664497) q[1];
sx q[1];
rz(-2.4058172) q[1];
sx q[1];
rz(1.9939569) q[1];
rz(-1.6699722) q[3];
sx q[3];
rz(-2.7226825) q[3];
sx q[3];
rz(0.3813972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2799985) q[2];
sx q[2];
rz(-1.7656606) q[2];
sx q[2];
rz(-0.17246788) q[2];
rz(2.5837768) q[3];
sx q[3];
rz(-2.6004801) q[3];
sx q[3];
rz(-1.7012168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43011343) q[0];
sx q[0];
rz(-2.5028296) q[0];
sx q[0];
rz(0.3399671) q[0];
rz(-2.4760447) q[1];
sx q[1];
rz(-1.7219209) q[1];
sx q[1];
rz(1.3131712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68379867) q[0];
sx q[0];
rz(-2.9556985) q[0];
sx q[0];
rz(-0.25769512) q[0];
rz(-pi) q[1];
rz(0.26076857) q[2];
sx q[2];
rz(-1.4968902) q[2];
sx q[2];
rz(0.36799612) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1024627) q[1];
sx q[1];
rz(-0.38929554) q[1];
sx q[1];
rz(-2.2401458) q[1];
x q[2];
rz(-1.4142939) q[3];
sx q[3];
rz(-1.1722574) q[3];
sx q[3];
rz(2.1074866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7790935) q[2];
sx q[2];
rz(-0.43467793) q[2];
sx q[2];
rz(-2.4707826) q[2];
rz(2.8166411) q[3];
sx q[3];
rz(-2.3281125) q[3];
sx q[3];
rz(-1.0130829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83298564) q[0];
sx q[0];
rz(-1.2815463) q[0];
sx q[0];
rz(1.4776342) q[0];
rz(-2.1113405) q[1];
sx q[1];
rz(-1.6696842) q[1];
sx q[1];
rz(-1.3546863) q[1];
rz(-0.75922913) q[2];
sx q[2];
rz(-0.6741796) q[2];
sx q[2];
rz(0.37311935) q[2];
rz(0.58086953) q[3];
sx q[3];
rz(-1.8884251) q[3];
sx q[3];
rz(-2.6407218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
