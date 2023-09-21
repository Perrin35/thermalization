OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.86413971) q[0];
sx q[0];
rz(-1.5530518) q[0];
sx q[0];
rz(-1.5074402) q[0];
rz(-1.545067) q[1];
sx q[1];
rz(-2.5453321) q[1];
sx q[1];
rz(2.526386) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9278487) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(2.5736789) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7841714) q[2];
sx q[2];
rz(-1.3961785) q[2];
sx q[2];
rz(-1.071196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33949172) q[1];
sx q[1];
rz(-1.1420982) q[1];
sx q[1];
rz(0.93356737) q[1];
x q[2];
rz(1.7513566) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(-0.80871049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7011828) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(-0.33828503) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(0.88589823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(0.030348226) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.679927) q[0];
sx q[0];
rz(-1.3711509) q[0];
sx q[0];
rz(-3.1398849) q[0];
rz(-pi) q[1];
x q[1];
rz(3.121071) q[2];
sx q[2];
rz(-1.1164718) q[2];
sx q[2];
rz(-0.12873912) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9873725) q[1];
sx q[1];
rz(-1.0013354) q[1];
sx q[1];
rz(-0.082055883) q[1];
rz(-pi) q[2];
rz(3.0475572) q[3];
sx q[3];
rz(-1.3841076) q[3];
sx q[3];
rz(-1.7969014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(-1.9937817) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(-0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(-2.8258064) q[0];
rz(2.2029927) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(-2.8895203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4323498) q[0];
sx q[0];
rz(-0.81626695) q[0];
sx q[0];
rz(-2.4791251) q[0];
x q[1];
rz(-2.4263072) q[2];
sx q[2];
rz(-0.89865696) q[2];
sx q[2];
rz(-0.31179024) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38110106) q[1];
sx q[1];
rz(-0.64844202) q[1];
sx q[1];
rz(1.8849424) q[1];
rz(-pi) q[2];
rz(-2.4113703) q[3];
sx q[3];
rz(-0.50958868) q[3];
sx q[3];
rz(-0.93917055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(0.034051731) q[2];
rz(3.1241336) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(-2.0461369) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(-1.6148286) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(-2.4345051) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93543816) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(0.48396707) q[0];
rz(-pi) q[1];
rz(2.884042) q[2];
sx q[2];
rz(-2.6817245) q[2];
sx q[2];
rz(2.3515153) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4428963) q[1];
sx q[1];
rz(-2.2203608) q[1];
sx q[1];
rz(0.15028468) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2624192) q[3];
sx q[3];
rz(-1.7763419) q[3];
sx q[3];
rz(-1.2984848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-2.6814931) q[2];
rz(1.7442616) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(2.0460515) q[1];
sx q[1];
rz(-1.6004326) q[1];
sx q[1];
rz(2.8869693) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1063227) q[0];
sx q[0];
rz(-1.584504) q[0];
sx q[0];
rz(-1.6499707) q[0];
rz(-pi) q[1];
rz(0.014572797) q[2];
sx q[2];
rz(-0.99327786) q[2];
sx q[2];
rz(0.84469675) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6756145) q[1];
sx q[1];
rz(-1.3295768) q[1];
sx q[1];
rz(0.66959186) q[1];
x q[2];
rz(-2.7084064) q[3];
sx q[3];
rz(-2.232589) q[3];
sx q[3];
rz(-1.2695241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.022481) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.508537) q[3];
sx q[3];
rz(-2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901533) q[0];
sx q[0];
rz(-2.4017161) q[0];
sx q[0];
rz(2.8421463) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(-2.9350231) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5305938) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(2.3321652) q[0];
x q[1];
rz(2.1913387) q[2];
sx q[2];
rz(-0.62180078) q[2];
sx q[2];
rz(-1.3602464) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5094604) q[1];
sx q[1];
rz(-1.7473979) q[1];
sx q[1];
rz(0.86061865) q[1];
rz(-pi) q[2];
rz(0.91512485) q[3];
sx q[3];
rz(-1.3279928) q[3];
sx q[3];
rz(-1.1806928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.841659) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(0.28277961) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(0.38152951) q[0];
rz(0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(1.8136224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8570003) q[0];
sx q[0];
rz(-0.45495957) q[0];
sx q[0];
rz(1.8332464) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46220025) q[2];
sx q[2];
rz(-2.1415347) q[2];
sx q[2];
rz(0.58909033) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5495587) q[1];
sx q[1];
rz(-0.70211071) q[1];
sx q[1];
rz(1.83105) q[1];
x q[2];
rz(-1.1442723) q[3];
sx q[3];
rz(-2.0445619) q[3];
sx q[3];
rz(1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61838377) q[2];
sx q[2];
rz(-2.0760459) q[2];
sx q[2];
rz(1.7810812) q[2];
rz(-1.4303738) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(-2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9946063) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(0.18187901) q[0];
rz(-0.47422844) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-0.95091933) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61688214) q[0];
sx q[0];
rz(-1.2843772) q[0];
sx q[0];
rz(-0.25047238) q[0];
rz(-pi) q[1];
rz(1.3457001) q[2];
sx q[2];
rz(-0.27268073) q[2];
sx q[2];
rz(2.8358012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7656895) q[1];
sx q[1];
rz(-1.6357058) q[1];
sx q[1];
rz(-2.9325571) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1675646) q[3];
sx q[3];
rz(-1.2207165) q[3];
sx q[3];
rz(-1.10266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9178847) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(0.63265911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52371812) q[0];
sx q[0];
rz(-1.0567559) q[0];
sx q[0];
rz(1.3762208) q[0];
rz(2.7245522) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(-0.65972796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8575681) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(2.2542623) q[0];
x q[1];
rz(-1.3525891) q[2];
sx q[2];
rz(-0.80112427) q[2];
sx q[2];
rz(1.4584691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8068741) q[1];
sx q[1];
rz(-1.8611307) q[1];
sx q[1];
rz(-0.90805407) q[1];
rz(-pi) q[2];
rz(-1.4755867) q[3];
sx q[3];
rz(-2.6009437) q[3];
sx q[3];
rz(-1.3117865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-1.0326577) q[3];
sx q[3];
rz(0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24348564) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.3056668) q[0];
rz(1.9650412) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(-1.0356888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0337692) q[0];
sx q[0];
rz(-2.6080837) q[0];
sx q[0];
rz(2.4643722) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56634855) q[2];
sx q[2];
rz(-0.95745917) q[2];
sx q[2];
rz(1.4373506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1199011) q[1];
sx q[1];
rz(-0.38858116) q[1];
sx q[1];
rz(-0.67967023) q[1];
rz(-pi) q[2];
rz(-0.74929897) q[3];
sx q[3];
rz(-1.370508) q[3];
sx q[3];
rz(-1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.62853652) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(-1.1516085) q[2];
rz(0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7619027) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(0.097461854) q[2];
sx q[2];
rz(-2.7290191) q[2];
sx q[2];
rz(1.4636427) q[2];
rz(1.5619754) q[3];
sx q[3];
rz(-2.1891441) q[3];
sx q[3];
rz(1.7563663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
