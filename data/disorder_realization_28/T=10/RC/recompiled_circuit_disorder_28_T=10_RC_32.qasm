OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2064535) q[0];
sx q[0];
rz(-0.78092617) q[0];
sx q[0];
rz(2.934802) q[0];
rz(-2.7517154) q[1];
sx q[1];
rz(-2.0808527) q[1];
sx q[1];
rz(-2.9614255) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8040112) q[0];
sx q[0];
rz(-2.0330226) q[0];
sx q[0];
rz(1.1929212) q[0];
x q[1];
rz(1.4397058) q[2];
sx q[2];
rz(-1.0616454) q[2];
sx q[2];
rz(-0.57123643) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85045056) q[1];
sx q[1];
rz(-2.2524815) q[1];
sx q[1];
rz(-3.046361) q[1];
x q[2];
rz(-1.7546685) q[3];
sx q[3];
rz(-2.2443218) q[3];
sx q[3];
rz(-2.9415188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41123286) q[2];
sx q[2];
rz(-0.87324548) q[2];
sx q[2];
rz(1.8475378) q[2];
rz(2.7358352) q[3];
sx q[3];
rz(-1.6399222) q[3];
sx q[3];
rz(-2.7348203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2186573) q[0];
sx q[0];
rz(-2.1056392) q[0];
sx q[0];
rz(-0.12582114) q[0];
rz(-0.80548349) q[1];
sx q[1];
rz(-0.80635726) q[1];
sx q[1];
rz(1.3719826) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94700891) q[0];
sx q[0];
rz(-2.2451375) q[0];
sx q[0];
rz(-1.6499004) q[0];
rz(0.54601045) q[2];
sx q[2];
rz(-1.464932) q[2];
sx q[2];
rz(0.0042303483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0224277) q[1];
sx q[1];
rz(-2.8875774) q[1];
sx q[1];
rz(-2.959842) q[1];
x q[2];
rz(-1.3268746) q[3];
sx q[3];
rz(-1.1141277) q[3];
sx q[3];
rz(-1.4671385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8877318) q[2];
sx q[2];
rz(-0.68887201) q[2];
sx q[2];
rz(-2.1453693) q[2];
rz(-2.0455202) q[3];
sx q[3];
rz(-1.6888065) q[3];
sx q[3];
rz(2.2912099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753733) q[0];
sx q[0];
rz(-2.1815364) q[0];
sx q[0];
rz(1.0590142) q[0];
rz(1.1478708) q[1];
sx q[1];
rz(-0.73906001) q[1];
sx q[1];
rz(2.0770729) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7721467) q[0];
sx q[0];
rz(-1.1755953) q[0];
sx q[0];
rz(-2.4482083) q[0];
rz(-pi) q[1];
rz(1.9469444) q[2];
sx q[2];
rz(-0.40824879) q[2];
sx q[2];
rz(0.86366913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73653446) q[1];
sx q[1];
rz(-1.5591803) q[1];
sx q[1];
rz(1.2606773) q[1];
x q[2];
rz(-0.16626658) q[3];
sx q[3];
rz(-2.2883752) q[3];
sx q[3];
rz(-2.1239514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0718096) q[2];
sx q[2];
rz(-1.9548543) q[2];
sx q[2];
rz(-0.26322571) q[2];
rz(1.1188544) q[3];
sx q[3];
rz(-1.2759195) q[3];
sx q[3];
rz(2.3222205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333106) q[0];
sx q[0];
rz(-0.044697035) q[0];
sx q[0];
rz(1.7472965) q[0];
rz(1.0385723) q[1];
sx q[1];
rz(-1.690052) q[1];
sx q[1];
rz(-2.0746453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79558668) q[0];
sx q[0];
rz(-2.6104402) q[0];
sx q[0];
rz(-1.4941494) q[0];
x q[1];
rz(-2.8550451) q[2];
sx q[2];
rz(-0.28038014) q[2];
sx q[2];
rz(-2.7718411) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9001138) q[1];
sx q[1];
rz(-0.80726868) q[1];
sx q[1];
rz(-0.57358731) q[1];
x q[2];
rz(2.9158343) q[3];
sx q[3];
rz(-2.0421713) q[3];
sx q[3];
rz(2.2798722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84714326) q[2];
sx q[2];
rz(-1.3501945) q[2];
sx q[2];
rz(-2.8520544) q[2];
rz(2.6211522) q[3];
sx q[3];
rz(-1.3402904) q[3];
sx q[3];
rz(-2.382544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4905869) q[0];
sx q[0];
rz(-3.1327972) q[0];
sx q[0];
rz(2.3068413) q[0];
rz(-1.897215) q[1];
sx q[1];
rz(-1.2507739) q[1];
sx q[1];
rz(1.429819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1236876) q[0];
sx q[0];
rz(-0.062739685) q[0];
sx q[0];
rz(0.040266589) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3313815) q[2];
sx q[2];
rz(-1.8481701) q[2];
sx q[2];
rz(1.1053567) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5809708) q[1];
sx q[1];
rz(-2.3498658) q[1];
sx q[1];
rz(-0.03246275) q[1];
rz(0.17794869) q[3];
sx q[3];
rz(-0.92015172) q[3];
sx q[3];
rz(-1.4072756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4804068) q[2];
sx q[2];
rz(-2.1358261) q[2];
sx q[2];
rz(1.8583813) q[2];
rz(-0.13218203) q[3];
sx q[3];
rz(-0.32871267) q[3];
sx q[3];
rz(0.33974084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64602393) q[0];
sx q[0];
rz(-0.060957242) q[0];
sx q[0];
rz(0.47750372) q[0];
rz(1.5006789) q[1];
sx q[1];
rz(-1.532282) q[1];
sx q[1];
rz(-0.57055155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43317023) q[0];
sx q[0];
rz(-0.47009429) q[0];
sx q[0];
rz(-0.8870468) q[0];
rz(-2.6887367) q[2];
sx q[2];
rz(-1.2315244) q[2];
sx q[2];
rz(-2.0138182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30299444) q[1];
sx q[1];
rz(-2.0717151) q[1];
sx q[1];
rz(2.3272728) q[1];
x q[2];
rz(-1.5930575) q[3];
sx q[3];
rz(-0.92261693) q[3];
sx q[3];
rz(-2.9472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70696124) q[2];
sx q[2];
rz(-1.0501477) q[2];
sx q[2];
rz(0.79664191) q[2];
rz(0.32026511) q[3];
sx q[3];
rz(-2.0551149) q[3];
sx q[3];
rz(-1.2789352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758078) q[0];
sx q[0];
rz(-1.0362352) q[0];
sx q[0];
rz(-0.35807034) q[0];
rz(0.2886731) q[1];
sx q[1];
rz(-0.52983785) q[1];
sx q[1];
rz(1.8310865) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0689439) q[0];
sx q[0];
rz(-1.7582969) q[0];
sx q[0];
rz(2.5995273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1475122) q[2];
sx q[2];
rz(-0.95816441) q[2];
sx q[2];
rz(-0.47555579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4431745) q[1];
sx q[1];
rz(-1.0647173) q[1];
sx q[1];
rz(1.0531823) q[1];
rz(-pi) q[2];
rz(0.86705039) q[3];
sx q[3];
rz(-0.82074814) q[3];
sx q[3];
rz(1.5914608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33401176) q[2];
sx q[2];
rz(-2.5568805) q[2];
sx q[2];
rz(0.99747783) q[2];
rz(2.5618662) q[3];
sx q[3];
rz(-2.4119191) q[3];
sx q[3];
rz(1.7060446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8758133) q[0];
sx q[0];
rz(-1.6515323) q[0];
sx q[0];
rz(2.8009801) q[0];
rz(1.9050725) q[1];
sx q[1];
rz(-2.3842936) q[1];
sx q[1];
rz(1.1901201) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2954138) q[0];
sx q[0];
rz(-0.23970397) q[0];
sx q[0];
rz(-1.833605) q[0];
rz(-pi) q[1];
rz(-0.068433381) q[2];
sx q[2];
rz(-1.6514196) q[2];
sx q[2];
rz(1.2863408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3548673) q[1];
sx q[1];
rz(-2.7495972) q[1];
sx q[1];
rz(0.74704945) q[1];
x q[2];
rz(0.90772273) q[3];
sx q[3];
rz(-2.472496) q[3];
sx q[3];
rz(0.95975403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.74165806) q[2];
sx q[2];
rz(-0.88230336) q[2];
sx q[2];
rz(-2.7271872) q[2];
rz(1.3828145) q[3];
sx q[3];
rz(-1.7410991) q[3];
sx q[3];
rz(0.32489052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8906616) q[0];
sx q[0];
rz(-1.0216167) q[0];
sx q[0];
rz(2.9898306) q[0];
rz(-1.7157308) q[1];
sx q[1];
rz(-1.3769923) q[1];
sx q[1];
rz(2.192416) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28578368) q[0];
sx q[0];
rz(-1.4306418) q[0];
sx q[0];
rz(-2.6967718) q[0];
x q[1];
rz(-1.7911712) q[2];
sx q[2];
rz(-0.95249635) q[2];
sx q[2];
rz(3.0439723) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4703428) q[1];
sx q[1];
rz(-1.8103231) q[1];
sx q[1];
rz(-2.4688979) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98032326) q[3];
sx q[3];
rz(-1.983641) q[3];
sx q[3];
rz(-2.1291898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40733797) q[2];
sx q[2];
rz(-2.7463425) q[2];
sx q[2];
rz(0.43241832) q[2];
rz(-1.6010823) q[3];
sx q[3];
rz(-1.6316905) q[3];
sx q[3];
rz(0.66175118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0703053) q[0];
sx q[0];
rz(-0.27305958) q[0];
sx q[0];
rz(0.37316698) q[0];
rz(-2.7846653) q[1];
sx q[1];
rz(-2.280805) q[1];
sx q[1];
rz(-2.4180791) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0821738) q[0];
sx q[0];
rz(-0.84616236) q[0];
sx q[0];
rz(-1.7436149) q[0];
rz(-pi) q[1];
rz(0.11156545) q[2];
sx q[2];
rz(-0.68250436) q[2];
sx q[2];
rz(-1.1139368) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.53612667) q[1];
sx q[1];
rz(-1.8408753) q[1];
sx q[1];
rz(-0.77401604) q[1];
rz(-pi) q[2];
rz(-2.8711653) q[3];
sx q[3];
rz(-1.575483) q[3];
sx q[3];
rz(-3.1386496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92131203) q[2];
sx q[2];
rz(-1.3833412) q[2];
sx q[2];
rz(2.5881361) q[2];
rz(0.71183318) q[3];
sx q[3];
rz(-2.7332941) q[3];
sx q[3];
rz(1.0420943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9217459) q[0];
sx q[0];
rz(-0.60315673) q[0];
sx q[0];
rz(-3.0043816) q[0];
rz(0.28221054) q[1];
sx q[1];
rz(-0.78411513) q[1];
sx q[1];
rz(0.93906739) q[1];
rz(-2.6454906) q[2];
sx q[2];
rz(-1.9719057) q[2];
sx q[2];
rz(3.0531648) q[2];
rz(2.5440352) q[3];
sx q[3];
rz(-0.47511027) q[3];
sx q[3];
rz(-2.6078754) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];