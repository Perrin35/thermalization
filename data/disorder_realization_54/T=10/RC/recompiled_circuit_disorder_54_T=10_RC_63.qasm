OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(2.8110992) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(5.4332241) q[1];
sx q[1];
rz(10.133893) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5753484) q[0];
sx q[0];
rz(-0.69501221) q[0];
sx q[0];
rz(-1.1027176) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4719226) q[2];
sx q[2];
rz(-0.31497248) q[2];
sx q[2];
rz(1.0613943) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.658537) q[1];
sx q[1];
rz(-2.4194948) q[1];
sx q[1];
rz(0.87941054) q[1];
x q[2];
rz(2.0907008) q[3];
sx q[3];
rz(-1.2399779) q[3];
sx q[3];
rz(-1.0370805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4101397) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(-0.62227917) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(0.91631779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8371671) q[0];
sx q[0];
rz(-2.2080748) q[0];
sx q[0];
rz(-2.5159555) q[0];
x q[1];
rz(-1.7158521) q[2];
sx q[2];
rz(-0.74160355) q[2];
sx q[2];
rz(-1.8650101) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.898145) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(-1.2136202) q[1];
rz(1.968859) q[3];
sx q[3];
rz(-1.5916087) q[3];
sx q[3];
rz(-1.2956937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-0.82143482) q[2];
rz(3.1243096) q[3];
sx q[3];
rz(-2.343943) q[3];
sx q[3];
rz(2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18773742) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(-1.0082555) q[0];
rz(3.1058274) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(0.52454138) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7099972) q[0];
sx q[0];
rz(-1.4072197) q[0];
sx q[0];
rz(3.1113935) q[0];
x q[1];
rz(-2.4630765) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(3.0845272) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9959065) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(-1.593959) q[1];
rz(-pi) q[2];
rz(1.9906304) q[3];
sx q[3];
rz(-1.435278) q[3];
sx q[3];
rz(0.47437048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77164578) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(1.6463722) q[2];
rz(-1.8203991) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-3.0526429) q[0];
rz(2.6308909) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(0.68960062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8423858) q[0];
sx q[0];
rz(-2.8028574) q[0];
sx q[0];
rz(2.6329106) q[0];
rz(-1.3454516) q[2];
sx q[2];
rz(-2.1225404) q[2];
sx q[2];
rz(-1.4622886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5213485) q[1];
sx q[1];
rz(-0.84940956) q[1];
sx q[1];
rz(2.0116624) q[1];
rz(-1.7161937) q[3];
sx q[3];
rz(-2.0510011) q[3];
sx q[3];
rz(0.2916382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66578635) q[2];
sx q[2];
rz(-2.4643354) q[2];
sx q[2];
rz(2.518667) q[2];
rz(1.1359435) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(0.583453) q[0];
rz(-1.1460229) q[1];
sx q[1];
rz(-1.6505516) q[1];
sx q[1];
rz(1.6437644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4125815) q[0];
sx q[0];
rz(-0.81775613) q[0];
sx q[0];
rz(-1.2110932) q[0];
rz(-pi) q[1];
rz(0.79767144) q[2];
sx q[2];
rz(-1.5016218) q[2];
sx q[2];
rz(0.56874146) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.076482) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-2.4379424) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51550128) q[3];
sx q[3];
rz(-1.4294129) q[3];
sx q[3];
rz(-3.0859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5636469) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(2.5031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435796) q[0];
sx q[0];
rz(-1.4587412) q[0];
sx q[0];
rz(3.0946099) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(1.7061589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2611321) q[0];
sx q[0];
rz(-1.4229703) q[0];
sx q[0];
rz(-0.79974215) q[0];
rz(-pi) q[1];
rz(-1.6666744) q[2];
sx q[2];
rz(-2.5898858) q[2];
sx q[2];
rz(0.72223896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3458851) q[1];
sx q[1];
rz(-2.7137623) q[1];
sx q[1];
rz(-1.2659628) q[1];
rz(-pi) q[2];
rz(-3.1141838) q[3];
sx q[3];
rz(-1.3331183) q[3];
sx q[3];
rz(2.4159367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-1.0152638) q[2];
sx q[2];
rz(3.0701239) q[2];
rz(-1.4525157) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(-0.92648363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(-1.8619934) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(1.0095899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6972835) q[0];
sx q[0];
rz(-0.59519207) q[0];
sx q[0];
rz(2.5213581) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7570653) q[2];
sx q[2];
rz(-1.3587917) q[2];
sx q[2];
rz(0.42696135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0040782) q[1];
sx q[1];
rz(-2.5234748) q[1];
sx q[1];
rz(2.532258) q[1];
x q[2];
rz(-2.5261577) q[3];
sx q[3];
rz(-0.86343599) q[3];
sx q[3];
rz(0.59188852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0499095) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(-0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(0.7827951) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(2.696864) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0895773) q[0];
sx q[0];
rz(-1.1451632) q[0];
sx q[0];
rz(-2.1968967) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3808001) q[2];
sx q[2];
rz(-2.3173601) q[2];
sx q[2];
rz(-2.9234718) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0900314) q[1];
sx q[1];
rz(-1.6981484) q[1];
sx q[1];
rz(-0.8831555) q[1];
x q[2];
rz(-2.8430311) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(-3.0446133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7198221) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(0.81531173) q[2];
rz(-1.0347962) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.3172654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8274882) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(2.8544193) q[0];
rz(2.9526967) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4591601) q[0];
sx q[0];
rz(-0.81572616) q[0];
sx q[0];
rz(1.1343603) q[0];
x q[1];
rz(3.095093) q[2];
sx q[2];
rz(-0.87247712) q[2];
sx q[2];
rz(-0.38682129) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13739535) q[1];
sx q[1];
rz(-1.9837556) q[1];
sx q[1];
rz(1.2910299) q[1];
rz(-pi) q[2];
rz(2.513047) q[3];
sx q[3];
rz(-2.5124031) q[3];
sx q[3];
rz(0.43720804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.93280783) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(-0.83958158) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.055450913) q[0];
sx q[0];
rz(-2.3738326) q[0];
sx q[0];
rz(1.0593876) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(0.25451452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8043038) q[0];
sx q[0];
rz(-0.24233195) q[0];
sx q[0];
rz(-1.2605577) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.326667) q[2];
sx q[2];
rz(-2.0460528) q[2];
sx q[2];
rz(2.313386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3855615) q[1];
sx q[1];
rz(-2.8061211) q[1];
sx q[1];
rz(0.24753333) q[1];
rz(-pi) q[2];
x q[2];
rz(1.47154) q[3];
sx q[3];
rz(-1.1343065) q[3];
sx q[3];
rz(2.1193159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(-2.0937031) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(-2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027325252) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(2.7813773) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(-2.9011177) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(-2.9110254) q[3];
sx q[3];
rz(-0.96944001) q[3];
sx q[3];
rz(2.159202) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
