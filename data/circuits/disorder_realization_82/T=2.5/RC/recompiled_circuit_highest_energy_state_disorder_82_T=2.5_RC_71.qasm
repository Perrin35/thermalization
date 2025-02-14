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
rz(3.7511711) q[0];
sx q[0];
rz(10.114976) q[0];
rz(-1.8742427) q[1];
sx q[1];
rz(6.7725362) q[1];
sx q[1];
rz(9.0803774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19669721) q[0];
sx q[0];
rz(-1.718194) q[0];
sx q[0];
rz(-1.7833254) q[0];
x q[1];
rz(2.6292388) q[2];
sx q[2];
rz(-2.6690581) q[2];
sx q[2];
rz(2.3290881) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59565249) q[1];
sx q[1];
rz(-1.3134688) q[1];
sx q[1];
rz(-0.20471065) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8841759) q[3];
sx q[3];
rz(-2.628577) q[3];
sx q[3];
rz(-3.0611924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.043618) q[2];
sx q[2];
rz(-2.2100885) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69720307) q[0];
sx q[0];
rz(-0.81962219) q[0];
sx q[0];
rz(-0.75253734) q[0];
rz(1.2731816) q[1];
sx q[1];
rz(-1.9877142) q[1];
sx q[1];
rz(0.85535991) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4580247) q[0];
sx q[0];
rz(-0.23580256) q[0];
sx q[0];
rz(1.1382415) q[0];
rz(-pi) q[1];
rz(-2.1925621) q[2];
sx q[2];
rz(-1.3675642) q[2];
sx q[2];
rz(-2.0469613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8458057) q[1];
sx q[1];
rz(-1.1870904) q[1];
sx q[1];
rz(1.7471261) q[1];
x q[2];
rz(3.0166808) q[3];
sx q[3];
rz(-0.20044573) q[3];
sx q[3];
rz(0.87504609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37318834) q[2];
sx q[2];
rz(-1.6104001) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7975467) q[0];
sx q[0];
rz(-2.5757289) q[0];
sx q[0];
rz(1.4669482) q[0];
rz(0.037874669) q[1];
sx q[1];
rz(-1.2118309) q[1];
sx q[1];
rz(-1.1143335) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9301611) q[0];
sx q[0];
rz(-0.4273779) q[0];
sx q[0];
rz(-1.2965802) q[0];
rz(-2.5882095) q[2];
sx q[2];
rz(-1.1107365) q[2];
sx q[2];
rz(-1.9529238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13521299) q[1];
sx q[1];
rz(-0.65875444) q[1];
sx q[1];
rz(-0.82243311) q[1];
x q[2];
rz(-1.1258026) q[3];
sx q[3];
rz(-2.5126406) q[3];
sx q[3];
rz(0.75643874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7797839) q[2];
sx q[2];
rz(-0.59541687) q[2];
sx q[2];
rz(-2.6105866) q[2];
rz(-0.94046721) q[3];
sx q[3];
rz(-0.85175025) q[3];
sx q[3];
rz(1.4550335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4130212) q[0];
sx q[0];
rz(-2.72609) q[0];
sx q[0];
rz(0.79237932) q[0];
rz(3.0089695) q[1];
sx q[1];
rz(-0.68999973) q[1];
sx q[1];
rz(2.2161868) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44418078) q[0];
sx q[0];
rz(-1.385293) q[0];
sx q[0];
rz(-2.9836487) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4508574) q[2];
sx q[2];
rz(-0.45443568) q[2];
sx q[2];
rz(-2.0023521) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.99605152) q[1];
sx q[1];
rz(-0.62623238) q[1];
sx q[1];
rz(1.224451) q[1];
rz(-pi) q[2];
rz(2.6230249) q[3];
sx q[3];
rz(-1.858497) q[3];
sx q[3];
rz(0.30687919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35599071) q[2];
sx q[2];
rz(-1.8428558) q[2];
sx q[2];
rz(-1.7388657) q[2];
rz(1.2518903) q[3];
sx q[3];
rz(-1.8391049) q[3];
sx q[3];
rz(-0.29786626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.70016015) q[0];
sx q[0];
rz(-2.1692363) q[0];
sx q[0];
rz(2.1527619) q[0];
rz(-1.5160457) q[1];
sx q[1];
rz(-1.4715618) q[1];
sx q[1];
rz(0.54738799) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1158041) q[0];
sx q[0];
rz(-2.0594547) q[0];
sx q[0];
rz(-1.5343094) q[0];
x q[1];
rz(0.36409521) q[2];
sx q[2];
rz(-1.4949706) q[2];
sx q[2];
rz(0.66616601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.84505397) q[1];
sx q[1];
rz(-2.2530139) q[1];
sx q[1];
rz(2.6107772) q[1];
rz(-pi) q[2];
rz(0.27606729) q[3];
sx q[3];
rz(-1.6652341) q[3];
sx q[3];
rz(-0.51968473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0872515) q[2];
sx q[2];
rz(-0.50830278) q[2];
sx q[2];
rz(2.9906719) q[2];
rz(-1.6357251) q[3];
sx q[3];
rz(-1.3003636) q[3];
sx q[3];
rz(-1.4668363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6557789) q[0];
sx q[0];
rz(-1.1968311) q[0];
sx q[0];
rz(2.6659513) q[0];
rz(0.42426839) q[1];
sx q[1];
rz(-1.905966) q[1];
sx q[1];
rz(-1.7686527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17775336) q[0];
sx q[0];
rz(-2.8256157) q[0];
sx q[0];
rz(0.58167268) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5504136) q[2];
sx q[2];
rz(-1.1496964) q[2];
sx q[2];
rz(-0.84122259) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13481278) q[1];
sx q[1];
rz(-2.105851) q[1];
sx q[1];
rz(-2.7538009) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8133393) q[3];
sx q[3];
rz(-2.6879592) q[3];
sx q[3];
rz(-2.2528474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0515392) q[2];
sx q[2];
rz(-0.97344437) q[2];
sx q[2];
rz(-3.0494087) q[2];
rz(-1.9696382) q[3];
sx q[3];
rz(-0.53297526) q[3];
sx q[3];
rz(2.2251825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.5343269) q[0];
sx q[0];
rz(-0.83616513) q[0];
sx q[0];
rz(1.0323866) q[0];
rz(3.0119925) q[1];
sx q[1];
rz(-1.7584636) q[1];
sx q[1];
rz(1.3547156) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.082026012) q[0];
sx q[0];
rz(-1.0903918) q[0];
sx q[0];
rz(-2.4405368) q[0];
rz(-2.3806055) q[2];
sx q[2];
rz(-0.47522173) q[2];
sx q[2];
rz(1.5409249) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9224841) q[1];
sx q[1];
rz(-2.0754083) q[1];
sx q[1];
rz(-2.497106) q[1];
x q[2];
rz(1.0951429) q[3];
sx q[3];
rz(-1.5009592) q[3];
sx q[3];
rz(1.8997518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20595343) q[2];
sx q[2];
rz(-2.8796068) q[2];
sx q[2];
rz(-1.6550753) q[2];
rz(2.4691811) q[3];
sx q[3];
rz(-1.0786062) q[3];
sx q[3];
rz(3.0730754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9392149) q[0];
sx q[0];
rz(-0.12513932) q[0];
sx q[0];
rz(2.8676721) q[0];
rz(-2.9778453) q[1];
sx q[1];
rz(-1.8293019) q[1];
sx q[1];
rz(-2.7408677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76213259) q[0];
sx q[0];
rz(-0.5665938) q[0];
sx q[0];
rz(-2.6521126) q[0];
x q[1];
rz(-1.2873307) q[2];
sx q[2];
rz(-2.7230883) q[2];
sx q[2];
rz(-2.037627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1355181) q[1];
sx q[1];
rz(-1.589314) q[1];
sx q[1];
rz(0.89864267) q[1];
rz(1.1464045) q[3];
sx q[3];
rz(-1.2514858) q[3];
sx q[3];
rz(0.13892787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87216941) q[2];
sx q[2];
rz(-0.17435208) q[2];
sx q[2];
rz(1.4738458) q[2];
rz(0.10146865) q[3];
sx q[3];
rz(-1.2227819) q[3];
sx q[3];
rz(1.7469223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2673016) q[0];
sx q[0];
rz(-2.3912781) q[0];
sx q[0];
rz(-0.0016203298) q[0];
rz(-2.5770309) q[1];
sx q[1];
rz(-1.8058585) q[1];
sx q[1];
rz(-0.50216215) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82715568) q[0];
sx q[0];
rz(-1.0183882) q[0];
sx q[0];
rz(2.3236426) q[0];
rz(-pi) q[1];
rz(2.5880626) q[2];
sx q[2];
rz(-1.2147673) q[2];
sx q[2];
rz(-1.219762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7495216) q[1];
sx q[1];
rz(-1.3765613) q[1];
sx q[1];
rz(0.51699395) q[1];
x q[2];
rz(1.2099482) q[3];
sx q[3];
rz(-0.54665414) q[3];
sx q[3];
rz(2.9846559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.042171176) q[2];
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
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5779483) q[0];
sx q[0];
rz(-0.36643323) q[0];
sx q[0];
rz(1.4165437) q[0];
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
rz(-1.6141629) q[0];
sx q[0];
rz(-0.50656453) q[0];
sx q[0];
rz(2.3247983) q[0];
rz(-pi) q[1];
rz(2.9014355) q[2];
sx q[2];
rz(-0.36061812) q[2];
sx q[2];
rz(-2.220253) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4135189) q[1];
sx q[1];
rz(-0.46183837) q[1];
sx q[1];
rz(-1.6175265) q[1];
x q[2];
rz(2.5429728) q[3];
sx q[3];
rz(-1.4685923) q[3];
sx q[3];
rz(-0.24598611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.766091) q[2];
sx q[2];
rz(-1.4644863) q[2];
sx q[2];
rz(-0.38140934) q[2];
rz(0.31960791) q[3];
sx q[3];
rz(-2.3242798) q[3];
sx q[3];
rz(-2.4735425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1565336) q[0];
sx q[0];
rz(-1.1548797) q[0];
sx q[0];
rz(-0.67697939) q[0];
rz(2.1398687) q[1];
sx q[1];
rz(-1.6698508) q[1];
sx q[1];
rz(2.225266) q[1];
rz(2.4763686) q[2];
sx q[2];
rz(-1.4736253) q[2];
sx q[2];
rz(1.8237154) q[2];
rz(-0.65883815) q[3];
sx q[3];
rz(-1.4821977) q[3];
sx q[3];
rz(2.341996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
