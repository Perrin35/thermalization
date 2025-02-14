OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8389799) q[0];
sx q[0];
rz(-1.7054727) q[0];
sx q[0];
rz(2.9247395) q[0];
rz(-0.4878374) q[1];
sx q[1];
rz(-0.87895972) q[1];
sx q[1];
rz(-0.15329696) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76591208) q[0];
sx q[0];
rz(-1.3049676) q[0];
sx q[0];
rz(1.2714809) q[0];
rz(2.3011123) q[2];
sx q[2];
rz(-2.4908713) q[2];
sx q[2];
rz(2.4706728) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5342321) q[1];
sx q[1];
rz(-2.5776754) q[1];
sx q[1];
rz(-0.47718559) q[1];
x q[2];
rz(-1.4073154) q[3];
sx q[3];
rz(-0.77115763) q[3];
sx q[3];
rz(-2.5389391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6903901) q[2];
sx q[2];
rz(-2.5610552) q[2];
sx q[2];
rz(0.85980493) q[2];
rz(0.23990038) q[3];
sx q[3];
rz(-2.4247215) q[3];
sx q[3];
rz(2.8587604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4848223) q[0];
sx q[0];
rz(-1.4161994) q[0];
sx q[0];
rz(-0.91944486) q[0];
rz(2.2942309) q[1];
sx q[1];
rz(-0.58440009) q[1];
sx q[1];
rz(-1.970361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7396512) q[0];
sx q[0];
rz(-1.8293132) q[0];
sx q[0];
rz(-0.052019428) q[0];
x q[1];
rz(2.809735) q[2];
sx q[2];
rz(-2.2584791) q[2];
sx q[2];
rz(0.0092384641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1619809) q[1];
sx q[1];
rz(-1.9165601) q[1];
sx q[1];
rz(-0.27372056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53816157) q[3];
sx q[3];
rz(-1.5892913) q[3];
sx q[3];
rz(-1.5539813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.072711572) q[2];
sx q[2];
rz(-2.2525807) q[2];
sx q[2];
rz(1.2129126) q[2];
rz(2.9291901) q[3];
sx q[3];
rz(-2.0272777) q[3];
sx q[3];
rz(-2.4685278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4251959) q[0];
sx q[0];
rz(-1.238751) q[0];
sx q[0];
rz(-0.58309251) q[0];
rz(1.8816226) q[1];
sx q[1];
rz(-0.59695736) q[1];
sx q[1];
rz(-1.8276385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010600032) q[0];
sx q[0];
rz(-2.6525462) q[0];
sx q[0];
rz(-2.2940192) q[0];
rz(2.407786) q[2];
sx q[2];
rz(-1.5467522) q[2];
sx q[2];
rz(-2.1847385) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3594234) q[1];
sx q[1];
rz(-0.12453989) q[1];
sx q[1];
rz(1.4119536) q[1];
rz(-1.9989955) q[3];
sx q[3];
rz(-1.7230265) q[3];
sx q[3];
rz(0.58109944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.49695697) q[2];
sx q[2];
rz(-0.76377112) q[2];
sx q[2];
rz(2.606707) q[2];
rz(0.48318091) q[3];
sx q[3];
rz(-1.6202241) q[3];
sx q[3];
rz(-3.0218637) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0907165) q[0];
sx q[0];
rz(-3.0829939) q[0];
sx q[0];
rz(1.6706049) q[0];
rz(-1.927467) q[1];
sx q[1];
rz(-1.7045226) q[1];
sx q[1];
rz(2.6953183) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31201115) q[0];
sx q[0];
rz(-2.952842) q[0];
sx q[0];
rz(-1.1971367) q[0];
rz(0.032164737) q[2];
sx q[2];
rz(-0.52903658) q[2];
sx q[2];
rz(0.47495237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91736356) q[1];
sx q[1];
rz(-2.434772) q[1];
sx q[1];
rz(-2.7429105) q[1];
x q[2];
rz(-2/(11*pi)) q[3];
sx q[3];
rz(-0.24316517) q[3];
sx q[3];
rz(-1.9609708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2176167) q[2];
sx q[2];
rz(-0.46102229) q[2];
sx q[2];
rz(0.32611845) q[2];
rz(0.41607949) q[3];
sx q[3];
rz(-1.5030428) q[3];
sx q[3];
rz(-0.7777586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.071082696) q[0];
sx q[0];
rz(-1.5776881) q[0];
sx q[0];
rz(-0.34410205) q[0];
rz(-2.7857419) q[1];
sx q[1];
rz(-0.75332037) q[1];
sx q[1];
rz(0.25486249) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5256115) q[0];
sx q[0];
rz(-0.92824844) q[0];
sx q[0];
rz(-1.2429953) q[0];
x q[1];
rz(-0.10019014) q[2];
sx q[2];
rz(-2.0311573) q[2];
sx q[2];
rz(3.0222297) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50824158) q[1];
sx q[1];
rz(-2.3741747) q[1];
sx q[1];
rz(-0.91752802) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9790824) q[3];
sx q[3];
rz(-0.84813839) q[3];
sx q[3];
rz(-1.4681787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7039589) q[2];
sx q[2];
rz(-2.2767229) q[2];
sx q[2];
rz(0.09058365) q[2];
rz(2.3352052) q[3];
sx q[3];
rz(-1.3017637) q[3];
sx q[3];
rz(-1.3069794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059747132) q[0];
sx q[0];
rz(-1.2223926) q[0];
sx q[0];
rz(-0.61862373) q[0];
rz(-2.5289358) q[1];
sx q[1];
rz(-0.63059348) q[1];
sx q[1];
rz(-2.9887066) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02462834) q[0];
sx q[0];
rz(-1.8061588) q[0];
sx q[0];
rz(-0.23925608) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15494187) q[2];
sx q[2];
rz(-0.58523387) q[2];
sx q[2];
rz(-2.3335148) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96079338) q[1];
sx q[1];
rz(-1.6241778) q[1];
sx q[1];
rz(1.175715) q[1];
rz(-0.024645667) q[3];
sx q[3];
rz(-2.1507015) q[3];
sx q[3];
rz(-0.49296311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9360518) q[2];
sx q[2];
rz(-0.80465332) q[2];
sx q[2];
rz(2.0118227) q[2];
rz(-1.0063082) q[3];
sx q[3];
rz(-1.3200503) q[3];
sx q[3];
rz(-1.4973076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048700843) q[0];
sx q[0];
rz(-1.0660271) q[0];
sx q[0];
rz(2.997828) q[0];
rz(2.9394506) q[1];
sx q[1];
rz(-0.527924) q[1];
sx q[1];
rz(1.082083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3976404) q[0];
sx q[0];
rz(-0.86371585) q[0];
sx q[0];
rz(2.6567064) q[0];
rz(-pi) q[1];
rz(0.73418243) q[2];
sx q[2];
rz(-0.95962822) q[2];
sx q[2];
rz(2.2847459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6134167) q[1];
sx q[1];
rz(-0.74494637) q[1];
sx q[1];
rz(1.8862392) q[1];
rz(-0.64855002) q[3];
sx q[3];
rz(-2.289384) q[3];
sx q[3];
rz(0.64283338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9988592) q[2];
sx q[2];
rz(-1.961668) q[2];
sx q[2];
rz(2.414523) q[2];
rz(-1.8266228) q[3];
sx q[3];
rz(-1.246779) q[3];
sx q[3];
rz(1.4962014) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8975248) q[0];
sx q[0];
rz(-2.0638564) q[0];
sx q[0];
rz(-1.9986073) q[0];
rz(0.22363981) q[1];
sx q[1];
rz(-0.63839212) q[1];
sx q[1];
rz(-1.9167831) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95668225) q[0];
sx q[0];
rz(-1.8326129) q[0];
sx q[0];
rz(-1.2774521) q[0];
rz(-pi) q[1];
rz(-0.46240004) q[2];
sx q[2];
rz(-0.77653304) q[2];
sx q[2];
rz(-2.3729216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8401007) q[1];
sx q[1];
rz(-0.21243851) q[1];
sx q[1];
rz(1.2919687) q[1];
rz(-pi) q[2];
rz(0.16815987) q[3];
sx q[3];
rz(-0.56795299) q[3];
sx q[3];
rz(-2.0310543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.58330047) q[2];
sx q[2];
rz(-1.6879098) q[2];
sx q[2];
rz(-0.80201403) q[2];
rz(-1.2174886) q[3];
sx q[3];
rz(-3.0018482) q[3];
sx q[3];
rz(1.8961934) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8121346) q[0];
sx q[0];
rz(-1.3690925) q[0];
sx q[0];
rz(1.580397) q[0];
rz(2.2691057) q[1];
sx q[1];
rz(-2.8552738) q[1];
sx q[1];
rz(2.7365541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0319388) q[0];
sx q[0];
rz(-2.4068916) q[0];
sx q[0];
rz(-0.74560179) q[0];
rz(-pi) q[1];
rz(-1.4680176) q[2];
sx q[2];
rz(-2.1609302) q[2];
sx q[2];
rz(2.0128606) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8499706) q[1];
sx q[1];
rz(-2.1287781) q[1];
sx q[1];
rz(-2.4172952) q[1];
x q[2];
rz(0.13073374) q[3];
sx q[3];
rz(-0.14821136) q[3];
sx q[3];
rz(-0.036202567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8672455) q[2];
sx q[2];
rz(-0.87928191) q[2];
sx q[2];
rz(2.530063) q[2];
rz(1.3036171) q[3];
sx q[3];
rz(-2.7795064) q[3];
sx q[3];
rz(-2.005579) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5030454) q[0];
sx q[0];
rz(-1.2940116) q[0];
sx q[0];
rz(0.053330388) q[0];
rz(-0.57180014) q[1];
sx q[1];
rz(-1.6245533) q[1];
sx q[1];
rz(-1.8624051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17634468) q[0];
sx q[0];
rz(-1.2851479) q[0];
sx q[0];
rz(1.3081416) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0903075) q[2];
sx q[2];
rz(-1.4391403) q[2];
sx q[2];
rz(-2.5327794) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1385344) q[1];
sx q[1];
rz(-1.6489677) q[1];
sx q[1];
rz(2.1581236) q[1];
rz(-pi) q[2];
rz(-1.199715) q[3];
sx q[3];
rz(-1.367192) q[3];
sx q[3];
rz(-0.18818391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.46644396) q[2];
sx q[2];
rz(-1.6089336) q[2];
sx q[2];
rz(2.4565728) q[2];
rz(0.78392616) q[3];
sx q[3];
rz(-2.7214366) q[3];
sx q[3];
rz(2.4043731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6017799) q[0];
sx q[0];
rz(-1.2417326) q[0];
sx q[0];
rz(-1.7057521) q[0];
rz(2.336179) q[1];
sx q[1];
rz(-0.46331159) q[1];
sx q[1];
rz(0.70809271) q[1];
rz(-2.5067301) q[2];
sx q[2];
rz(-2.786953) q[2];
sx q[2];
rz(2.0155596) q[2];
rz(1.4516674) q[3];
sx q[3];
rz(-1.1584251) q[3];
sx q[3];
rz(2.3746015) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
