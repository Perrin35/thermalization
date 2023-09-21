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
rz(2.3072825) q[0];
sx q[0];
rz(6.351525) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(0.037820427) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8605182) q[0];
sx q[0];
rz(-1.3776508) q[0];
sx q[0];
rz(3.0618969) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90934609) q[2];
sx q[2];
rz(-1.6644018) q[2];
sx q[2];
rz(-1.0974761) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6112069) q[1];
sx q[1];
rz(-1.621765) q[1];
sx q[1];
rz(-0.74688046) q[1];
rz(-pi) q[2];
rz(-1.171265) q[3];
sx q[3];
rz(-2.7683308) q[3];
sx q[3];
rz(-1.8217306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35090703) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(-1.8581871) q[2];
rz(0.12617271) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(0.090601966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.693817) q[0];
sx q[0];
rz(-2.2651146) q[0];
sx q[0];
rz(2.5449261) q[0];
rz(1.5860575) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(1.7780875) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3421905) q[0];
sx q[0];
rz(-1.0418833) q[0];
sx q[0];
rz(2.0687813) q[0];
rz(-pi) q[1];
rz(-0.081625799) q[2];
sx q[2];
rz(-0.55742369) q[2];
sx q[2];
rz(0.73478414) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1494257) q[1];
sx q[1];
rz(-2.4352695) q[1];
sx q[1];
rz(-2.7041433) q[1];
rz(2.4035461) q[3];
sx q[3];
rz(-0.9934721) q[3];
sx q[3];
rz(0.20418024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3229225) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(-0.73454109) q[2];
rz(2.5143886) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(-2.396615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7221786) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(2.869379) q[0];
rz(-2.294337) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(-0.31262696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0756829) q[0];
sx q[0];
rz(-1.624375) q[0];
sx q[0];
rz(-2.9220394) q[0];
rz(-1.1001415) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(-0.91980308) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0215175) q[1];
sx q[1];
rz(-2.8719963) q[1];
sx q[1];
rz(1.3070379) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3787141) q[3];
sx q[3];
rz(-1.8490013) q[3];
sx q[3];
rz(-1.1080081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12038885) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(-2.3266501) q[2];
rz(-0.96873823) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(-0.43280861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75354904) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(-0.62336212) q[0];
rz(0.81758824) q[1];
sx q[1];
rz(-2.8249884) q[1];
sx q[1];
rz(-3.0923016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8636297) q[0];
sx q[0];
rz(-2.3290714) q[0];
sx q[0];
rz(0.066594007) q[0];
rz(-0.31844278) q[2];
sx q[2];
rz(-2.2568984) q[2];
sx q[2];
rz(-2.7903914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66344391) q[1];
sx q[1];
rz(-1.899292) q[1];
sx q[1];
rz(2.8847787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6922675) q[3];
sx q[3];
rz(-2.0530564) q[3];
sx q[3];
rz(1.1376017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4108882) q[2];
sx q[2];
rz(-1.5894019) q[2];
sx q[2];
rz(-2.1566186) q[2];
rz(2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(2.8682017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7966998) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(3.0539736) q[0];
rz(-0.15631974) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(2.2706251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6773416) q[0];
sx q[0];
rz(-0.28872492) q[0];
sx q[0];
rz(0.10084734) q[0];
x q[1];
rz(1.2147374) q[2];
sx q[2];
rz(-0.41741727) q[2];
sx q[2];
rz(2.2603214) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37967967) q[1];
sx q[1];
rz(-1.8508136) q[1];
sx q[1];
rz(0.67358394) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17467588) q[3];
sx q[3];
rz(-2.3274765) q[3];
sx q[3];
rz(2.0811618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.093420371) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(2.2757163) q[3];
sx q[3];
rz(-2.4398949) q[3];
sx q[3];
rz(0.29233366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6784994) q[0];
sx q[0];
rz(-2.0149639) q[0];
sx q[0];
rz(-3.034806) q[0];
rz(-1.9550025) q[1];
sx q[1];
rz(-2.1148966) q[1];
sx q[1];
rz(-1.0245163) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5750676) q[0];
sx q[0];
rz(-1.1199513) q[0];
sx q[0];
rz(-1.2498115) q[0];
x q[1];
rz(2.5249135) q[2];
sx q[2];
rz(-2.4384535) q[2];
sx q[2];
rz(-0.48381915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5974821) q[1];
sx q[1];
rz(-2.1253617) q[1];
sx q[1];
rz(-2.5764562) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5173353) q[3];
sx q[3];
rz(-0.58371021) q[3];
sx q[3];
rz(2.7492085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(0.8708896) q[2];
rz(-1.332256) q[3];
sx q[3];
rz(-0.94227666) q[3];
sx q[3];
rz(1.4107305) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1720599) q[0];
sx q[0];
rz(-1.4071858) q[0];
sx q[0];
rz(-2.5906738) q[0];
rz(2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(-0.23682061) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0222586) q[0];
sx q[0];
rz(-3.0354781) q[0];
sx q[0];
rz(-1.9342213) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4332438) q[2];
sx q[2];
rz(-1.7195065) q[2];
sx q[2];
rz(-0.48352048) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6432861) q[1];
sx q[1];
rz(-2.8611538) q[1];
sx q[1];
rz(1.2966869) q[1];
rz(-1.2348433) q[3];
sx q[3];
rz(-0.7630322) q[3];
sx q[3];
rz(-1.2515765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44234309) q[2];
sx q[2];
rz(-2.4599059) q[2];
sx q[2];
rz(-0.44011763) q[2];
rz(-2.0464499) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(1.346689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31297627) q[0];
sx q[0];
rz(-1.799311) q[0];
sx q[0];
rz(-3.112088) q[0];
rz(-0.94738952) q[1];
sx q[1];
rz(-0.82059971) q[1];
sx q[1];
rz(-3.0227919) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5598181) q[0];
sx q[0];
rz(-1.8811328) q[0];
sx q[0];
rz(-1.52542) q[0];
rz(-pi) q[1];
x q[1];
rz(2.115909) q[2];
sx q[2];
rz(-1.9695749) q[2];
sx q[2];
rz(1.1023956) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6606632) q[1];
sx q[1];
rz(-2.5813563) q[1];
sx q[1];
rz(2.5653097) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16468594) q[3];
sx q[3];
rz(-1.0822923) q[3];
sx q[3];
rz(0.96795852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3999195) q[2];
sx q[2];
rz(-2.6779149) q[2];
sx q[2];
rz(-1.5681533) q[2];
rz(-1.9741156) q[3];
sx q[3];
rz(-1.4195331) q[3];
sx q[3];
rz(-1.8673816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90299273) q[0];
sx q[0];
rz(-2.3165343) q[0];
sx q[0];
rz(-0.15326823) q[0];
rz(-1.0614456) q[1];
sx q[1];
rz(-2.5669284) q[1];
sx q[1];
rz(1.1538039) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35837367) q[0];
sx q[0];
rz(-1.8405387) q[0];
sx q[0];
rz(-2.7620035) q[0];
rz(-0.88887631) q[2];
sx q[2];
rz(-2.0098915) q[2];
sx q[2];
rz(-3.0438434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.02281636) q[1];
sx q[1];
rz(-2.9252242) q[1];
sx q[1];
rz(-1.8323891) q[1];
rz(-pi) q[2];
rz(3.0204569) q[3];
sx q[3];
rz(-2.2545358) q[3];
sx q[3];
rz(-2.5072806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8447421) q[2];
sx q[2];
rz(-1.1721609) q[2];
sx q[2];
rz(1.6142169) q[2];
rz(2.1447003) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(2.2629288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24755724) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(0.19009185) q[0];
rz(-2.4709573) q[1];
sx q[1];
rz(-1.9452483) q[1];
sx q[1];
rz(1.6533096) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9207536) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(0.12551813) q[0];
rz(-pi) q[1];
rz(1.3004296) q[2];
sx q[2];
rz(-1.9989982) q[2];
sx q[2];
rz(2.2141475) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8390159) q[1];
sx q[1];
rz(-0.95514983) q[1];
sx q[1];
rz(2.2018196) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58191396) q[3];
sx q[3];
rz(-0.96393185) q[3];
sx q[3];
rz(-2.5556263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0037447475) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(-0.89912644) q[2];
rz(2.6265465) q[3];
sx q[3];
rz(-2.6656272) q[3];
sx q[3];
rz(-0.17764828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0777733) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(2.8521815) q[1];
sx q[1];
rz(-2.2402973) q[1];
sx q[1];
rz(-1.4351861) q[1];
rz(-2.8171956) q[2];
sx q[2];
rz(-0.73642052) q[2];
sx q[2];
rz(-3.0036075) q[2];
rz(1.8437456) q[3];
sx q[3];
rz(-1.6932586) q[3];
sx q[3];
rz(0.3441588) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
