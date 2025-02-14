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
rz(-0.43871969) q[0];
sx q[0];
rz(-2.5320142) q[0];
sx q[0];
rz(0.69019812) q[0];
rz(-1.8742427) q[1];
sx q[1];
rz(-2.6522418) q[1];
sx q[1];
rz(-2.7971921) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448954) q[0];
sx q[0];
rz(-1.4233986) q[0];
sx q[0];
rz(1.3582673) q[0];
rz(-2.6292388) q[2];
sx q[2];
rz(-0.47253451) q[2];
sx q[2];
rz(-0.81250459) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8616383) q[1];
sx q[1];
rz(-0.32740232) q[1];
sx q[1];
rz(2.2285217) q[1];
rz(2.9696483) q[3];
sx q[3];
rz(-2.056582) q[3];
sx q[3];
rz(-2.7051089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0979746) q[2];
sx q[2];
rz(-0.93150413) q[2];
sx q[2];
rz(-2.0841058) q[2];
rz(-0.99772325) q[3];
sx q[3];
rz(-1.7505373) q[3];
sx q[3];
rz(-1.9690431) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69720307) q[0];
sx q[0];
rz(-0.81962219) q[0];
sx q[0];
rz(-2.3890553) q[0];
rz(1.8684111) q[1];
sx q[1];
rz(-1.9877142) q[1];
sx q[1];
rz(2.2862327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1269) q[0];
sx q[0];
rz(-1.3570667) q[0];
sx q[0];
rz(0.10038169) q[0];
rz(-pi) q[1];
rz(-0.94903058) q[2];
sx q[2];
rz(-1.7740284) q[2];
sx q[2];
rz(-2.0469613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8458057) q[1];
sx q[1];
rz(-1.9545022) q[1];
sx q[1];
rz(-1.3944666) q[1];
rz(-2.9426676) q[3];
sx q[3];
rz(-1.595605) q[3];
sx q[3];
rz(2.5682784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7684043) q[2];
sx q[2];
rz(-1.5311925) q[2];
sx q[2];
rz(-1.9739523) q[2];
rz(3.043637) q[3];
sx q[3];
rz(-0.28147134) q[3];
sx q[3];
rz(-1.9907985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(-0.34404594) q[0];
sx q[0];
rz(-2.5757289) q[0];
sx q[0];
rz(1.6746445) q[0];
rz(3.103718) q[1];
sx q[1];
rz(-1.9297618) q[1];
sx q[1];
rz(2.0272592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0328465) q[0];
sx q[0];
rz(-1.4583197) q[0];
sx q[0];
rz(-1.9839909) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5882095) q[2];
sx q[2];
rz(-2.0308562) q[2];
sx q[2];
rz(-1.1886688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73015651) q[1];
sx q[1];
rz(-1.1056285) q[1];
sx q[1];
rz(-0.48480715) q[1];
rz(-2.15184) q[3];
sx q[3];
rz(-1.3147745) q[3];
sx q[3];
rz(2.6953002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36180878) q[2];
sx q[2];
rz(-2.5461758) q[2];
sx q[2];
rz(2.6105866) q[2];
rz(2.2011254) q[3];
sx q[3];
rz(-0.85175025) q[3];
sx q[3];
rz(1.4550335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4130212) q[0];
sx q[0];
rz(-2.72609) q[0];
sx q[0];
rz(-2.3492133) q[0];
rz(-3.0089695) q[1];
sx q[1];
rz(-2.4515929) q[1];
sx q[1];
rz(-0.92540583) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6974119) q[0];
sx q[0];
rz(-1.385293) q[0];
sx q[0];
rz(-2.9836487) q[0];
rz(1.4508574) q[2];
sx q[2];
rz(-2.687157) q[2];
sx q[2];
rz(-2.0023521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99605152) q[1];
sx q[1];
rz(-2.5153603) q[1];
sx q[1];
rz(-1.9171417) q[1];
x q[2];
rz(-2.6230249) q[3];
sx q[3];
rz(-1.2830956) q[3];
sx q[3];
rz(0.30687919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7856019) q[2];
sx q[2];
rz(-1.8428558) q[2];
sx q[2];
rz(1.7388657) q[2];
rz(-1.8897024) q[3];
sx q[3];
rz(-1.8391049) q[3];
sx q[3];
rz(2.8437264) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70016015) q[0];
sx q[0];
rz(-2.1692363) q[0];
sx q[0];
rz(2.1527619) q[0];
rz(1.6255469) q[1];
sx q[1];
rz(-1.6700309) q[1];
sx q[1];
rz(-0.54738799) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0382017) q[0];
sx q[0];
rz(-2.6516838) q[0];
sx q[0];
rz(0.068504917) q[0];
rz(-2.9314032) q[2];
sx q[2];
rz(-2.7700305) q[2];
sx q[2];
rz(0.7084058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.095371072) q[1];
sx q[1];
rz(-2.3041572) q[1];
sx q[1];
rz(1.0135021) q[1];
x q[2];
rz(-1.6689273) q[3];
sx q[3];
rz(-1.2959912) q[3];
sx q[3];
rz(2.0637728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0872515) q[2];
sx q[2];
rz(-2.6332899) q[2];
sx q[2];
rz(-0.15092078) q[2];
rz(-1.5058676) q[3];
sx q[3];
rz(-1.8412291) q[3];
sx q[3];
rz(1.6747564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48581377) q[0];
sx q[0];
rz(-1.1968311) q[0];
sx q[0];
rz(-2.6659513) q[0];
rz(-2.7173243) q[1];
sx q[1];
rz(-1.2356267) q[1];
sx q[1];
rz(1.7686527) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1899446) q[0];
sx q[0];
rz(-1.7423672) q[0];
sx q[0];
rz(-2.8749332) q[0];
rz(-pi) q[1];
rz(0.67694725) q[2];
sx q[2];
rz(-2.430655) q[2];
sx q[2];
rz(-2.9589911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9109505) q[1];
sx q[1];
rz(-1.2394718) q[1];
sx q[1];
rz(2.1403007) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0128355) q[3];
sx q[3];
rz(-1.6762432) q[3];
sx q[3];
rz(2.6783708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0515392) q[2];
sx q[2];
rz(-0.97344437) q[2];
sx q[2];
rz(-3.0494087) q[2];
rz(1.1719545) q[3];
sx q[3];
rz(-0.53297526) q[3];
sx q[3];
rz(-0.91641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6072657) q[0];
sx q[0];
rz(-0.83616513) q[0];
sx q[0];
rz(-2.109206) q[0];
rz(3.0119925) q[1];
sx q[1];
rz(-1.383129) q[1];
sx q[1];
rz(1.786877) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.860703) q[0];
sx q[0];
rz(-0.96183813) q[0];
sx q[0];
rz(2.1693098) q[0];
x q[1];
rz(2.7849136) q[2];
sx q[2];
rz(-1.891815) q[2];
sx q[2];
rz(-0.73250801) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70020218) q[1];
sx q[1];
rz(-1.017015) q[1];
sx q[1];
rz(-0.9662083) q[1];
rz(0.078523648) q[3];
sx q[3];
rz(-2.0451945) q[3];
sx q[3];
rz(-0.36488381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.20595343) q[2];
sx q[2];
rz(-0.26198584) q[2];
sx q[2];
rz(1.6550753) q[2];
rz(-2.4691811) q[3];
sx q[3];
rz(-2.0629864) q[3];
sx q[3];
rz(3.0730754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9392149) q[0];
sx q[0];
rz(-3.0164533) q[0];
sx q[0];
rz(-2.8676721) q[0];
rz(-2.9778453) q[1];
sx q[1];
rz(-1.3122908) q[1];
sx q[1];
rz(-0.40072498) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3253097) q[0];
sx q[0];
rz(-2.0643215) q[0];
sx q[0];
rz(1.8614344) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.854262) q[2];
sx q[2];
rz(-2.7230883) q[2];
sx q[2];
rz(-1.1039656) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0060746) q[1];
sx q[1];
rz(-1.5522787) q[1];
sx q[1];
rz(2.24295) q[1];
rz(-pi) q[2];
rz(2.2473281) q[3];
sx q[3];
rz(-0.52519631) q[3];
sx q[3];
rz(2.0390455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87216941) q[2];
sx q[2];
rz(-2.9672406) q[2];
sx q[2];
rz(-1.4738458) q[2];
rz(3.040124) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(-1.7469223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87429109) q[0];
sx q[0];
rz(-2.3912781) q[0];
sx q[0];
rz(3.1399723) q[0];
rz(0.56456176) q[1];
sx q[1];
rz(-1.3357342) q[1];
sx q[1];
rz(-2.6394305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9084887) q[0];
sx q[0];
rz(-0.90051578) q[0];
sx q[0];
rz(-2.3044808) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55353005) q[2];
sx q[2];
rz(-1.2147673) q[2];
sx q[2];
rz(-1.219762) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6356023) q[1];
sx q[1];
rz(-0.5491623) q[1];
sx q[1];
rz(-2.7628187) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21162866) q[3];
sx q[3];
rz(-1.0628848) q[3];
sx q[3];
rz(-2.882559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0994215) q[2];
sx q[2];
rz(-1.9440938) q[2];
sx q[2];
rz(-1.8966804) q[2];
rz(0.20243195) q[3];
sx q[3];
rz(-1.3916241) q[3];
sx q[3];
rz(2.1421053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56364432) q[0];
sx q[0];
rz(-2.7751594) q[0];
sx q[0];
rz(1.7250489) q[0];
rz(-1.1997403) q[1];
sx q[1];
rz(-1.9681294) q[1];
sx q[1];
rz(2.8660668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107127) q[0];
sx q[0];
rz(-1.2322324) q[0];
sx q[0];
rz(1.186446) q[0];
rz(-1.4813384) q[2];
sx q[2];
rz(-1.2209792) q[2];
sx q[2];
rz(-1.9642771) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72807377) q[1];
sx q[1];
rz(-2.6797543) q[1];
sx q[1];
rz(1.6175265) q[1];
rz(2.5429728) q[3];
sx q[3];
rz(-1.4685923) q[3];
sx q[3];
rz(-0.24598611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3755017) q[2];
sx q[2];
rz(-1.6771064) q[2];
sx q[2];
rz(-0.38140934) q[2];
rz(-2.8219847) q[3];
sx q[3];
rz(-2.3242798) q[3];
sx q[3];
rz(-2.4735425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9850591) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(-1.001724) q[1];
sx q[1];
rz(-1.6698508) q[1];
sx q[1];
rz(2.225266) q[1];
rz(-2.4763686) q[2];
sx q[2];
rz(-1.6679674) q[2];
sx q[2];
rz(-1.3178772) q[2];
rz(2.4827545) q[3];
sx q[3];
rz(-1.4821977) q[3];
sx q[3];
rz(2.341996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
