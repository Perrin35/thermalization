OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(4.1783279) q[1];
sx q[1];
rz(8.7611603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27204313) q[0];
sx q[0];
rz(-2.6798195) q[0];
sx q[0];
rz(3.0884577) q[0];
rz(-0.77779777) q[2];
sx q[2];
rz(-1.8282837) q[2];
sx q[2];
rz(0.68629005) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44890468) q[1];
sx q[1];
rz(-1.9069888) q[1];
sx q[1];
rz(0.2775788) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1575559) q[3];
sx q[3];
rz(-0.26502702) q[3];
sx q[3];
rz(0.36039263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-3.0905241) q[2];
rz(2.5845394) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(-1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58650815) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(0.82582981) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85428836) q[0];
sx q[0];
rz(-2.0137631) q[0];
sx q[0];
rz(1.7914377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77483564) q[2];
sx q[2];
rz(-0.92445395) q[2];
sx q[2];
rz(-1.5409842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.021124161) q[1];
sx q[1];
rz(-1.3509343) q[1];
sx q[1];
rz(-3.1411509) q[1];
x q[2];
rz(1.5561043) q[3];
sx q[3];
rz(-2.5689295) q[3];
sx q[3];
rz(-0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.4276918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817292) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(0.088009134) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2217076) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(2.1973781) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4918409) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(0.31313716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.50476915) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.4298646) q[1];
rz(-pi) q[2];
rz(3.0089278) q[3];
sx q[3];
rz(-2.3972315) q[3];
sx q[3];
rz(-2.6877407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(0.50764817) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(-1.5159336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.078331) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(-2.9177833) q[0];
rz(1.0117202) q[2];
sx q[2];
rz(-1.498073) q[2];
sx q[2];
rz(0.59414547) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8856636) q[1];
sx q[1];
rz(-2.4293578) q[1];
sx q[1];
rz(0.42339143) q[1];
rz(2.5528615) q[3];
sx q[3];
rz(-1.7541459) q[3];
sx q[3];
rz(2.6453032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97840727) q[0];
sx q[0];
rz(-0.95046959) q[0];
sx q[0];
rz(1.0582256) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72563719) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(-0.63647905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27069651) q[1];
sx q[1];
rz(-2.0082698) q[1];
sx q[1];
rz(-2.2938964) q[1];
rz(-pi) q[2];
rz(-0.026167913) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(-0.22931306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58364761) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(-1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(0.46498743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2842448) q[0];
sx q[0];
rz(-1.2372969) q[0];
sx q[0];
rz(-2.0302982) q[0];
x q[1];
rz(-0.34157413) q[2];
sx q[2];
rz(-1.6401059) q[2];
sx q[2];
rz(-0.14788936) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2496693) q[1];
sx q[1];
rz(-1.5230471) q[1];
sx q[1];
rz(2.2062917) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24352169) q[3];
sx q[3];
rz(-2.4176819) q[3];
sx q[3];
rz(-0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(0.93368357) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(-0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(1.4200462) q[0];
rz(3.1177915) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(-2.9856317) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1056846) q[0];
sx q[0];
rz(-1.7419635) q[0];
sx q[0];
rz(0.25733421) q[0];
rz(-pi) q[1];
rz(-0.51257001) q[2];
sx q[2];
rz(-1.9480431) q[2];
sx q[2];
rz(2.9104428) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.29282001) q[1];
sx q[1];
rz(-2.9584868) q[1];
sx q[1];
rz(-1.8715026) q[1];
x q[2];
rz(-0.74219269) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.069313958) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(1.0726661) q[2];
rz(-0.32564751) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(0.33777133) q[0];
rz(-2.0514964) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-0.24857323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45390689) q[0];
sx q[0];
rz(-2.0810063) q[0];
sx q[0];
rz(1.8574255) q[0];
rz(-2.3209004) q[2];
sx q[2];
rz(-0.65260115) q[2];
sx q[2];
rz(1.5479659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(-2.2680125) q[1];
rz(-3.0646938) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(-0.83991915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(-3.0548813) q[2];
rz(0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865737) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-2.8588262) q[0];
rz(-0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(-1.3185906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0193664) q[0];
sx q[0];
rz(-2.0113693) q[0];
sx q[0];
rz(0.088935436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65638541) q[2];
sx q[2];
rz(-2.0771386) q[2];
sx q[2];
rz(-2.2040747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88398891) q[1];
sx q[1];
rz(-2.249243) q[1];
sx q[1];
rz(0.18737327) q[1];
rz(-pi) q[2];
rz(-0.75026476) q[3];
sx q[3];
rz(-1.8599469) q[3];
sx q[3];
rz(0.17444785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7302154) q[2];
sx q[2];
rz(-2.8432379) q[2];
sx q[2];
rz(0.39917699) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(-1.5426853) q[0];
rz(-2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.261196) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0364089) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(1.4559792) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5485498) q[2];
sx q[2];
rz(-0.97264475) q[2];
sx q[2];
rz(-1.0215789) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0295769) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(-3.0081248) q[1];
rz(-pi) q[2];
rz(-0.19461467) q[3];
sx q[3];
rz(-2.448304) q[3];
sx q[3];
rz(-0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5579055) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-0.56979257) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(2.5789554) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(-1.530004) q[2];
sx q[2];
rz(-0.57208021) q[2];
sx q[2];
rz(-0.013442599) q[2];
rz(-2.5504997) q[3];
sx q[3];
rz(-1.686284) q[3];
sx q[3];
rz(-1.4808663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
