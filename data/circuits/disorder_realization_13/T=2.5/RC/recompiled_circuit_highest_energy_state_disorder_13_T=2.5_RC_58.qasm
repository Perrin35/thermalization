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
rz(-2.8615992) q[0];
sx q[0];
rz(-2.1092829) q[0];
sx q[0];
rz(-1.9814459) q[0];
rz(2.2634444) q[1];
sx q[1];
rz(-0.4385837) q[1];
sx q[1];
rz(-2.261472) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11186799) q[0];
sx q[0];
rz(-2.1639185) q[0];
sx q[0];
rz(-0.2082227) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7574667) q[2];
sx q[2];
rz(-2.7323033) q[2];
sx q[2];
rz(-0.14044204) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82950006) q[1];
sx q[1];
rz(-2.1662427) q[1];
sx q[1];
rz(-0.19705806) q[1];
rz(-2.9008242) q[3];
sx q[3];
rz(-2.6759522) q[3];
sx q[3];
rz(-2.9064399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0264414) q[2];
sx q[2];
rz(-2.5796964) q[2];
sx q[2];
rz(2.0217516) q[2];
rz(1.0037054) q[3];
sx q[3];
rz(-0.34990889) q[3];
sx q[3];
rz(0.85163918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.411946) q[0];
sx q[0];
rz(-0.8929407) q[0];
sx q[0];
rz(-2.6935284) q[0];
rz(-1.0753151) q[1];
sx q[1];
rz(-2.0649464) q[1];
sx q[1];
rz(-3.101128) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69175289) q[0];
sx q[0];
rz(-1.3794071) q[0];
sx q[0];
rz(-2.2474849) q[0];
rz(-1.1574817) q[2];
sx q[2];
rz(-1.4885474) q[2];
sx q[2];
rz(2.8678107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76772122) q[1];
sx q[1];
rz(-1.5295856) q[1];
sx q[1];
rz(2.8727864) q[1];
rz(-pi) q[2];
rz(-2.8289757) q[3];
sx q[3];
rz(-1.5626057) q[3];
sx q[3];
rz(0.20144486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.86541092) q[2];
sx q[2];
rz(-2.2409596) q[2];
sx q[2];
rz(2.6382228) q[2];
rz(-1.6940176) q[3];
sx q[3];
rz(-1.2332799) q[3];
sx q[3];
rz(-2.2342822) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68509787) q[0];
sx q[0];
rz(-0.95360294) q[0];
sx q[0];
rz(-2.8431235) q[0];
rz(1.5553156) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(1.5063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28373805) q[0];
sx q[0];
rz(-1.2728134) q[0];
sx q[0];
rz(-1.4922754) q[0];
x q[1];
rz(-1.9048018) q[2];
sx q[2];
rz(-2.4270505) q[2];
sx q[2];
rz(1.3151907) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8019234) q[1];
sx q[1];
rz(-2.443772) q[1];
sx q[1];
rz(2.9342447) q[1];
rz(-pi) q[2];
rz(1.8465145) q[3];
sx q[3];
rz(-0.19681588) q[3];
sx q[3];
rz(-2.4829602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35206587) q[2];
sx q[2];
rz(-2.3730998) q[2];
sx q[2];
rz(-2.3699769) q[2];
rz(1.4113034) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(2.3119149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4378433) q[0];
sx q[0];
rz(-0.65173906) q[0];
sx q[0];
rz(1.334345) q[0];
rz(1.6288545) q[1];
sx q[1];
rz(-1.6213497) q[1];
sx q[1];
rz(-0.13660647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2481987) q[0];
sx q[0];
rz(-1.2635487) q[0];
sx q[0];
rz(0.91715468) q[0];
rz(-1.8804276) q[2];
sx q[2];
rz(-1.6694434) q[2];
sx q[2];
rz(-1.6473687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1699511) q[1];
sx q[1];
rz(-1.882269) q[1];
sx q[1];
rz(-0.71572742) q[1];
rz(-pi) q[2];
rz(-2.3312373) q[3];
sx q[3];
rz(-2.6763066) q[3];
sx q[3];
rz(-2.3734059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85285774) q[2];
sx q[2];
rz(-1.4790269) q[2];
sx q[2];
rz(-2.9020201) q[2];
rz(2.1824956) q[3];
sx q[3];
rz(-1.164091) q[3];
sx q[3];
rz(-2.1036928) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94944823) q[0];
sx q[0];
rz(-1.0553772) q[0];
sx q[0];
rz(1.5978093) q[0];
rz(1.1100618) q[1];
sx q[1];
rz(-1.203457) q[1];
sx q[1];
rz(-1.69453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2758165) q[0];
sx q[0];
rz(-1.5779243) q[0];
sx q[0];
rz(2.843186) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1008312) q[2];
sx q[2];
rz(-1.9136179) q[2];
sx q[2];
rz(-1.1072323) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6441783) q[1];
sx q[1];
rz(-1.7458054) q[1];
sx q[1];
rz(-0.67703663) q[1];
x q[2];
rz(2.5980159) q[3];
sx q[3];
rz(-1.2050306) q[3];
sx q[3];
rz(0.32777946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1453104) q[2];
sx q[2];
rz(-2.0029533) q[2];
sx q[2];
rz(-0.81926695) q[2];
rz(2.1938358) q[3];
sx q[3];
rz(-0.78815931) q[3];
sx q[3];
rz(1.2558827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.509023) q[0];
sx q[0];
rz(-0.13795723) q[0];
sx q[0];
rz(-2.5914958) q[0];
rz(2.8944648) q[1];
sx q[1];
rz(-1.988966) q[1];
sx q[1];
rz(0.94608847) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9086015) q[0];
sx q[0];
rz(-2.6938022) q[0];
sx q[0];
rz(-0.93771387) q[0];
rz(-pi) q[1];
rz(-1.7839512) q[2];
sx q[2];
rz(-1.6317131) q[2];
sx q[2];
rz(-2.3062381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4934686) q[1];
sx q[1];
rz(-1.5350684) q[1];
sx q[1];
rz(2.6717527) q[1];
x q[2];
rz(-0.40762679) q[3];
sx q[3];
rz(-2.8823311) q[3];
sx q[3];
rz(-1.2490937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3863824) q[2];
sx q[2];
rz(-0.55042616) q[2];
sx q[2];
rz(0.68106252) q[2];
rz(-0.81124535) q[3];
sx q[3];
rz(-1.0439876) q[3];
sx q[3];
rz(-2.5913141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1193467) q[0];
sx q[0];
rz(-0.62048727) q[0];
sx q[0];
rz(-2.3959809) q[0];
rz(-1.9781808) q[1];
sx q[1];
rz(-2.5201576) q[1];
sx q[1];
rz(0.17722873) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772607) q[0];
sx q[0];
rz(-1.2998616) q[0];
sx q[0];
rz(1.9146108) q[0];
rz(-pi) q[1];
rz(0.96259768) q[2];
sx q[2];
rz(-0.11044914) q[2];
sx q[2];
rz(-0.77753528) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9504576) q[1];
sx q[1];
rz(-2.0759575) q[1];
sx q[1];
rz(2.7209366) q[1];
rz(0.050595119) q[3];
sx q[3];
rz(-0.6953763) q[3];
sx q[3];
rz(-0.68948244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54997286) q[2];
sx q[2];
rz(-0.98229304) q[2];
sx q[2];
rz(-1.4493235) q[2];
rz(-2.5189404) q[3];
sx q[3];
rz(-0.52716523) q[3];
sx q[3];
rz(-1.8221633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016668884) q[0];
sx q[0];
rz(-1.7542087) q[0];
sx q[0];
rz(1.5119875) q[0];
rz(0.92388693) q[1];
sx q[1];
rz(-1.7216564) q[1];
sx q[1];
rz(0.61663827) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409326) q[0];
sx q[0];
rz(-0.97447878) q[0];
sx q[0];
rz(2.9572972) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.086313649) q[2];
sx q[2];
rz(-0.46937916) q[2];
sx q[2];
rz(-2.1522107) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5376327) q[1];
sx q[1];
rz(-0.62509552) q[1];
sx q[1];
rz(2.9949466) q[1];
x q[2];
rz(1.8455335) q[3];
sx q[3];
rz(-1.0062075) q[3];
sx q[3];
rz(0.72455154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7743249) q[2];
sx q[2];
rz(-0.35949817) q[2];
sx q[2];
rz(0.16505879) q[2];
rz(0.74602357) q[3];
sx q[3];
rz(-1.6227928) q[3];
sx q[3];
rz(0.042796854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54166722) q[0];
sx q[0];
rz(-2.8494819) q[0];
sx q[0];
rz(2.7274729) q[0];
rz(-2.7060624) q[1];
sx q[1];
rz(-0.51594096) q[1];
sx q[1];
rz(-1.2783277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1762487) q[0];
sx q[0];
rz(-0.2834191) q[0];
sx q[0];
rz(-2.118628) q[0];
rz(-pi) q[1];
rz(-0.31164557) q[2];
sx q[2];
rz(-1.7362744) q[2];
sx q[2];
rz(-0.32725016) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2486677) q[1];
sx q[1];
rz(-1.5589412) q[1];
sx q[1];
rz(3.0737033) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0986365) q[3];
sx q[3];
rz(-1.5164638) q[3];
sx q[3];
rz(-2.899037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4182959) q[2];
sx q[2];
rz(-1.3822684) q[2];
sx q[2];
rz(0.74271512) q[2];
rz(0.79139477) q[3];
sx q[3];
rz(-2.4580038) q[3];
sx q[3];
rz(1.7323823) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0940229) q[0];
sx q[0];
rz(-2.8122734) q[0];
sx q[0];
rz(2.2299715) q[0];
rz(-0.98512638) q[1];
sx q[1];
rz(-0.42675012) q[1];
sx q[1];
rz(1.8034579) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3365173) q[0];
sx q[0];
rz(-2.0455845) q[0];
sx q[0];
rz(2.7874448) q[0];
x q[1];
rz(-2.2608032) q[2];
sx q[2];
rz(-2.6652209) q[2];
sx q[2];
rz(2.7169189) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9175901) q[1];
sx q[1];
rz(-2.3907967) q[1];
sx q[1];
rz(-2.177813) q[1];
rz(-pi) q[2];
rz(0.33439016) q[3];
sx q[3];
rz(-1.2161126) q[3];
sx q[3];
rz(-2.3324764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7868339) q[2];
sx q[2];
rz(-1.4749196) q[2];
sx q[2];
rz(-1.9762529) q[2];
rz(1.745584) q[3];
sx q[3];
rz(-1.952012) q[3];
sx q[3];
rz(1.5497807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-1.6373445) q[0];
sx q[0];
rz(-1.5865542) q[0];
sx q[0];
rz(-2.4028548) q[0];
rz(-0.63738102) q[1];
sx q[1];
rz(-1.5370054) q[1];
sx q[1];
rz(-0.76269033) q[1];
rz(-1.9945831) q[2];
sx q[2];
rz(-1.4374566) q[2];
sx q[2];
rz(-0.79604436) q[2];
rz(-1.1203587) q[3];
sx q[3];
rz(-1.0584581) q[3];
sx q[3];
rz(-1.2367482) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
