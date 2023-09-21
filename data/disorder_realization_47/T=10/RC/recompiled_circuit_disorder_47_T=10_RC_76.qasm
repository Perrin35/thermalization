OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(3.0732529) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(0.037820427) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8671994) q[0];
sx q[0];
rz(-1.4925856) q[0];
sx q[0];
rz(1.377051) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0231789) q[2];
sx q[2];
rz(-2.2288433) q[2];
sx q[2];
rz(2.7409035) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99327982) q[1];
sx q[1];
rz(-0.82511745) q[1];
sx q[1];
rz(-1.5013904) q[1];
rz(-pi) q[2];
rz(-1.9703276) q[3];
sx q[3];
rz(-0.37326187) q[3];
sx q[3];
rz(1.319862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7906856) q[2];
sx q[2];
rz(-2.6280792) q[2];
sx q[2];
rz(-1.8581871) q[2];
rz(0.12617271) q[3];
sx q[3];
rz(-1.7254555) q[3];
sx q[3];
rz(-0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44777563) q[0];
sx q[0];
rz(-2.2651146) q[0];
sx q[0];
rz(2.5449261) q[0];
rz(-1.5860575) q[1];
sx q[1];
rz(-1.7998453) q[1];
sx q[1];
rz(1.7780875) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024433) q[0];
sx q[0];
rz(-1.1457232) q[0];
sx q[0];
rz(-0.58702472) q[0];
rz(-pi) q[1];
rz(-2.5856651) q[2];
sx q[2];
rz(-1.613942) q[2];
sx q[2];
rz(2.3748929) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9047782) q[1];
sx q[1];
rz(-1.8493376) q[1];
sx q[1];
rz(0.65794557) q[1];
rz(-0.76916839) q[3];
sx q[3];
rz(-0.90220074) q[3];
sx q[3];
rz(1.234496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8186701) q[2];
sx q[2];
rz(-1.0007891) q[2];
sx q[2];
rz(0.73454109) q[2];
rz(-0.62720403) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(-2.396615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(-0.27221361) q[0];
rz(-0.84725562) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(-2.8289657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659098) q[0];
sx q[0];
rz(-1.5172177) q[0];
sx q[0];
rz(-2.9220394) q[0];
rz(-pi) q[1];
rz(-2.59358) q[2];
sx q[2];
rz(-1.9805038) q[2];
sx q[2];
rz(0.41248413) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9885013) q[1];
sx q[1];
rz(-1.830849) q[1];
sx q[1];
rz(-0.071916332) q[1];
rz(2.3787141) q[3];
sx q[3];
rz(-1.8490013) q[3];
sx q[3];
rz(1.1080081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0212038) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(-2.1728544) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75354904) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(-0.62336212) q[0];
rz(-0.81758824) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-0.049291074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76694861) q[0];
sx q[0];
rz(-2.3809803) q[0];
sx q[0];
rz(1.6409372) q[0];
rz(0.31844278) q[2];
sx q[2];
rz(-0.88469425) q[2];
sx q[2];
rz(-2.7903914) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.019776736) q[1];
sx q[1];
rz(-2.7275118) q[1];
sx q[1];
rz(0.93036923) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44932511) q[3];
sx q[3];
rz(-2.0530564) q[3];
sx q[3];
rz(-2.0039909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4108882) q[2];
sx q[2];
rz(-1.5521908) q[2];
sx q[2];
rz(-0.98497406) q[2];
rz(-0.90304053) q[3];
sx q[3];
rz(-1.9786381) q[3];
sx q[3];
rz(0.27339098) q[3];
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
rz(-2.7966998) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(-0.087619089) q[0];
rz(0.15631974) q[1];
sx q[1];
rz(-0.55115288) q[1];
sx q[1];
rz(2.2706251) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5721711) q[0];
sx q[0];
rz(-1.2835802) q[0];
sx q[0];
rz(-1.6006908) q[0];
rz(-1.1768612) q[2];
sx q[2];
rz(-1.7125868) q[2];
sx q[2];
rz(-2.1243492) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2838382) q[1];
sx q[1];
rz(-0.72099599) q[1];
sx q[1];
rz(-2.7096219) q[1];
x q[2];
rz(-1.3887651) q[3];
sx q[3];
rz(-2.3689299) q[3];
sx q[3];
rz(-2.3327737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.093420371) q[2];
sx q[2];
rz(-0.8478567) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(-2.2757163) q[3];
sx q[3];
rz(-0.70169774) q[3];
sx q[3];
rz(-2.849259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784994) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(-0.10678664) q[0];
rz(-1.1865901) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(-1.0245163) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5750676) q[0];
sx q[0];
rz(-1.1199513) q[0];
sx q[0];
rz(1.2498115) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60501955) q[2];
sx q[2];
rz(-1.1875249) q[2];
sx q[2];
rz(-0.59125102) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7926327) q[1];
sx q[1];
rz(-2.0434725) q[1];
sx q[1];
rz(2.2036168) q[1];
x q[2];
rz(-0.6242574) q[3];
sx q[3];
rz(-0.58371021) q[3];
sx q[3];
rz(-2.7492085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7496877) q[2];
sx q[2];
rz(-1.9275894) q[2];
sx q[2];
rz(2.270703) q[2];
rz(-1.8093367) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(1.7308621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(2.5906738) q[0];
rz(-0.46539601) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-0.23682061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.119334) q[0];
sx q[0];
rz(-0.10611457) q[0];
sx q[0];
rz(1.2073713) q[0];
rz(1.7083488) q[2];
sx q[2];
rz(-1.7195065) q[2];
sx q[2];
rz(-2.6580722) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6432861) q[1];
sx q[1];
rz(-2.8611538) q[1];
sx q[1];
rz(-1.8449057) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83644609) q[3];
sx q[3];
rz(-1.340938) q[3];
sx q[3];
rz(-3.0695855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6992496) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(2.701475) q[2];
rz(-2.0464499) q[3];
sx q[3];
rz(-1.4705855) q[3];
sx q[3];
rz(1.7949036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8286164) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(-3.112088) q[0];
rz(-2.1942031) q[1];
sx q[1];
rz(-2.3209929) q[1];
sx q[1];
rz(0.11880076) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122075) q[0];
sx q[0];
rz(-2.8280624) q[0];
sx q[0];
rz(-0.14051147) q[0];
rz(-pi) q[1];
rz(0.45784874) q[2];
sx q[2];
rz(-2.0689788) q[2];
sx q[2];
rz(0.23715167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6606632) q[1];
sx q[1];
rz(-0.56023635) q[1];
sx q[1];
rz(0.57628298) q[1];
rz(1.8700065) q[3];
sx q[3];
rz(-0.51338235) q[3];
sx q[3];
rz(-1.8332924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3999195) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(1.5734394) q[2];
rz(1.167477) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(0.15326823) q[0];
rz(1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.9877888) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3183347) q[0];
sx q[0];
rz(-1.2055921) q[0];
sx q[0];
rz(1.860114) q[0];
x q[1];
rz(2.2527163) q[2];
sx q[2];
rz(-2.0098915) q[2];
sx q[2];
rz(-3.0438434) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.24473083) q[1];
sx q[1];
rz(-1.7796928) q[1];
sx q[1];
rz(0.056785866) q[1];
x q[2];
rz(-1.7180175) q[3];
sx q[3];
rz(-2.4489093) q[3];
sx q[3];
rz(2.6976531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.29685059) q[2];
sx q[2];
rz(-1.9694318) q[2];
sx q[2];
rz(1.5273757) q[2];
rz(0.99689233) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24755724) q[0];
sx q[0];
rz(-2.0443125) q[0];
sx q[0];
rz(-0.19009185) q[0];
rz(0.67063531) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(-1.488283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3538441) q[0];
sx q[0];
rz(-1.4453381) q[0];
sx q[0];
rz(1.6018484) q[0];
rz(0.52942099) q[2];
sx q[2];
rz(-2.6396857) q[2];
sx q[2];
rz(1.6255962) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8390159) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(-0.93977309) q[1];
rz(-2.2641803) q[3];
sx q[3];
rz(-1.1023695) q[3];
sx q[3];
rz(1.7978158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1378479) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(-2.2424662) q[2];
rz(-0.51504618) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(0.17764828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0638194) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(-2.8521815) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(-0.70984798) q[2];
sx q[2];
rz(-1.786543) q[2];
sx q[2];
rz(-1.67698) q[2];
rz(0.12712052) q[3];
sx q[3];
rz(-1.2999429) q[3];
sx q[3];
rz(-1.1924549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
