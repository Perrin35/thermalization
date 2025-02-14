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
rz(0.76437104) q[0];
sx q[0];
rz(-1.3410913) q[0];
sx q[0];
rz(-0.90547019) q[0];
rz(4.6484923) q[1];
sx q[1];
rz(5.3223106) q[1];
sx q[1];
rz(7.9237908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95661058) q[0];
sx q[0];
rz(-0.8763823) q[0];
sx q[0];
rz(2.882769) q[0];
rz(-0.17619422) q[2];
sx q[2];
rz(-1.1374047) q[2];
sx q[2];
rz(-0.36461634) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5792306) q[1];
sx q[1];
rz(-0.048431245) q[1];
sx q[1];
rz(-2.6391451) q[1];
x q[2];
rz(-2.1330058) q[3];
sx q[3];
rz(-1.8671745) q[3];
sx q[3];
rz(2.5135771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2291439) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(-0.32192117) q[2];
rz(-0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(1.9979075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8973273) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(-2.6296997) q[0];
rz(1.5882675) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(-2.0194676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.666709) q[0];
sx q[0];
rz(-0.32829744) q[0];
sx q[0];
rz(2.8216381) q[0];
rz(-2.000359) q[2];
sx q[2];
rz(-1.6074174) q[2];
sx q[2];
rz(-3.1319194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1074688) q[1];
sx q[1];
rz(-1.144088) q[1];
sx q[1];
rz(0.36860768) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6038325) q[3];
sx q[3];
rz(-2.3590292) q[3];
sx q[3];
rz(-1.7478706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7684218) q[2];
sx q[2];
rz(-1.827652) q[2];
sx q[2];
rz(-2.6785417) q[2];
rz(-2.566346) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(0.099253207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333882) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(-0.43011618) q[0];
rz(-0.46562132) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(0.63708416) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4642536) q[0];
sx q[0];
rz(-1.5298944) q[0];
sx q[0];
rz(0.014057191) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3561959) q[2];
sx q[2];
rz(-1.5893418) q[2];
sx q[2];
rz(-0.77814276) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3780669) q[1];
sx q[1];
rz(-0.77516205) q[1];
sx q[1];
rz(-2.7833392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1807084) q[3];
sx q[3];
rz(-2.0566787) q[3];
sx q[3];
rz(2.7220243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.78274) q[2];
sx q[2];
rz(-2.095486) q[2];
sx q[2];
rz(-0.72506881) q[2];
rz(-1.3310165) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(-0.63841188) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8722039) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(1.4440906) q[0];
rz(-0.048642453) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(2.8526502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.671593) q[0];
sx q[0];
rz(-0.44737383) q[0];
sx q[0];
rz(2.1567221) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3861395) q[2];
sx q[2];
rz(-0.15310213) q[2];
sx q[2];
rz(-0.36400041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.98499456) q[1];
sx q[1];
rz(-2.9356476) q[1];
sx q[1];
rz(-2.1316281) q[1];
rz(-pi) q[2];
rz(-0.66529243) q[3];
sx q[3];
rz(-1.6136618) q[3];
sx q[3];
rz(2.8957518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.284953) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(-2.7613769) q[2];
rz(-0.63816655) q[3];
sx q[3];
rz(-1.9117982) q[3];
sx q[3];
rz(2.0130472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3087092) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(-2.047245) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(-2.0424776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35661429) q[0];
sx q[0];
rz(-1.7081424) q[0];
sx q[0];
rz(2.7524968) q[0];
rz(-pi) q[1];
rz(2.3540456) q[2];
sx q[2];
rz(-0.468245) q[2];
sx q[2];
rz(-2.5494273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.079541072) q[1];
sx q[1];
rz(-1.7289484) q[1];
sx q[1];
rz(2.1248716) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3640251) q[3];
sx q[3];
rz(-0.66477699) q[3];
sx q[3];
rz(-0.95486508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8984453) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(0.97664991) q[2];
rz(-2.6099033) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(2.4272052) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0831182) q[0];
sx q[0];
rz(-2.0396905) q[0];
sx q[0];
rz(1.8699159) q[0];
rz(-2.7187128) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(1.5333102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9031808) q[0];
sx q[0];
rz(-0.067771284) q[0];
sx q[0];
rz(1.3890024) q[0];
x q[1];
rz(1.4752059) q[2];
sx q[2];
rz(-0.098473452) q[2];
sx q[2];
rz(-2.5661039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.21580885) q[1];
sx q[1];
rz(-1.6263576) q[1];
sx q[1];
rz(2.8381398) q[1];
rz(-pi) q[2];
rz(2.3284376) q[3];
sx q[3];
rz(-1.7405563) q[3];
sx q[3];
rz(-0.90857279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5292458) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(2.7563654) q[2];
rz(0.084608229) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(-2.2658074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2200634) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(-2.7440942) q[0];
rz(0.53783224) q[1];
sx q[1];
rz(-2.718524) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5631112) q[0];
sx q[0];
rz(-2.3627776) q[0];
sx q[0];
rz(-1.7527197) q[0];
rz(2.053431) q[2];
sx q[2];
rz(-0.74649278) q[2];
sx q[2];
rz(1.6255524) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7518378) q[1];
sx q[1];
rz(-1.0026649) q[1];
sx q[1];
rz(-2.0271432) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3621906) q[3];
sx q[3];
rz(-0.85242295) q[3];
sx q[3];
rz(1.8441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6323382) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(0.21305591) q[2];
rz(-2.3325855) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(1.2808778) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0179366) q[0];
sx q[0];
rz(-1.2024711) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(-2.842438) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(-2.1655653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3487745) q[0];
sx q[0];
rz(-2.8716757) q[0];
sx q[0];
rz(-2.7151373) q[0];
rz(2.9891122) q[2];
sx q[2];
rz(-1.7465542) q[2];
sx q[2];
rz(2.3703142) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.9240771) q[1];
sx q[1];
rz(-1.2189606) q[1];
sx q[1];
rz(0.25434964) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7175432) q[3];
sx q[3];
rz(-1.3433787) q[3];
sx q[3];
rz(0.018674803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8396847) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(-2.8470993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(-1.4991722) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(1.3444208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62468707) q[0];
sx q[0];
rz(-0.4962099) q[0];
sx q[0];
rz(-0.93006706) q[0];
rz(-0.70602472) q[2];
sx q[2];
rz(-2.386552) q[2];
sx q[2];
rz(-1.6598827) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0186179) q[1];
sx q[1];
rz(-2.5135165) q[1];
sx q[1];
rz(2.0043892) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2559861) q[3];
sx q[3];
rz(-0.84784283) q[3];
sx q[3];
rz(-0.9566488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53238791) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(-1.6877635) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5880599) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.4916627) q[0];
rz(-2.5904169) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(-2.5490882) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1051286) q[0];
sx q[0];
rz(-1.291664) q[0];
sx q[0];
rz(-1.3238086) q[0];
rz(-pi) q[1];
rz(0.76982381) q[2];
sx q[2];
rz(-0.69902674) q[2];
sx q[2];
rz(-1.5848643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3687835) q[1];
sx q[1];
rz(-1.9516489) q[1];
sx q[1];
rz(-1.8541149) q[1];
x q[2];
rz(0.09162147) q[3];
sx q[3];
rz(-0.50096782) q[3];
sx q[3];
rz(-2.9660781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7117915) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(0.05833021) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534828) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(1.7755605) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(0.57353061) q[2];
sx q[2];
rz(-1.7885359) q[2];
sx q[2];
rz(0.89319695) q[2];
rz(2.1128863) q[3];
sx q[3];
rz(-2.0301314) q[3];
sx q[3];
rz(0.83740656) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
