OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.70541731) q[0];
sx q[0];
rz(-2.5751312) q[0];
sx q[0];
rz(2.9705272) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(1.8106102) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0199466) q[0];
sx q[0];
rz(-2.2197147) q[0];
sx q[0];
rz(-1.5050423) q[0];
x q[1];
rz(0.51585977) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(3.0992103) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.060974412) q[1];
sx q[1];
rz(-1.3197101) q[1];
sx q[1];
rz(0.017647839) q[1];
x q[2];
rz(0.29403789) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(-0.77120632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(1.8120871) q[2];
rz(-1.7261516) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(0.20234385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(1.6888899) q[0];
rz(-2.9406722) q[1];
sx q[1];
rz(-1.0849489) q[1];
sx q[1];
rz(-0.30028775) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072104134) q[0];
sx q[0];
rz(-2.527643) q[0];
sx q[0];
rz(2.9314562) q[0];
rz(-pi) q[1];
rz(1.8882206) q[2];
sx q[2];
rz(-1.023264) q[2];
sx q[2];
rz(2.4878793) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2806432) q[1];
sx q[1];
rz(-1.7935392) q[1];
sx q[1];
rz(1.2692979) q[1];
rz(-pi) q[2];
rz(-2.977936) q[3];
sx q[3];
rz(-1.7573342) q[3];
sx q[3];
rz(1.2638307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(-0.84233061) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(2.8252576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91352275) q[0];
sx q[0];
rz(-1.5726349) q[0];
sx q[0];
rz(1.3257922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5306273) q[2];
sx q[2];
rz(-2.4749711) q[2];
sx q[2];
rz(1.1952458) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9005147) q[1];
sx q[1];
rz(-1.9652275) q[1];
sx q[1];
rz(-0.48836744) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29052492) q[3];
sx q[3];
rz(-1.2447262) q[3];
sx q[3];
rz(2.3323004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1091653) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(-0.12604776) q[3];
sx q[3];
rz(-1.3879317) q[3];
sx q[3];
rz(1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(-0.8316935) q[0];
rz(-1.7968934) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.5136738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154402) q[0];
sx q[0];
rz(-0.8988131) q[0];
sx q[0];
rz(-1.5120904) q[0];
x q[1];
rz(1.6985745) q[2];
sx q[2];
rz(-2.5737692) q[2];
sx q[2];
rz(-0.77726269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1308243) q[1];
sx q[1];
rz(-1.1225892) q[1];
sx q[1];
rz(0.018647714) q[1];
x q[2];
rz(-1.1553331) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(2.8727939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4251129) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(-0.41695693) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(2.6127889) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083374627) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(1.8378687) q[0];
rz(-0.45571348) q[1];
sx q[1];
rz(-0.2625176) q[1];
sx q[1];
rz(-0.051503332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542959) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(-3.0412578) q[0];
rz(1.8565882) q[2];
sx q[2];
rz(-1.4225866) q[2];
sx q[2];
rz(-2.1697901) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81899535) q[1];
sx q[1];
rz(-1.7261191) q[1];
sx q[1];
rz(1.1197234) q[1];
rz(-pi) q[2];
rz(-0.13637654) q[3];
sx q[3];
rz(-0.73835056) q[3];
sx q[3];
rz(-2.7969558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.7752825) q[2];
sx q[2];
rz(-0.59801897) q[2];
rz(2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0710058) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(-2.9396074) q[0];
rz(2.1760991) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-0.083267033) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88887963) q[0];
sx q[0];
rz(-1.3759334) q[0];
sx q[0];
rz(-1.3937508) q[0];
x q[1];
rz(-2.8667738) q[2];
sx q[2];
rz(-1.8289369) q[2];
sx q[2];
rz(-1.2061662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80379936) q[1];
sx q[1];
rz(-2.7593136) q[1];
sx q[1];
rz(-1.8610958) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6827624) q[3];
sx q[3];
rz(-0.60744951) q[3];
sx q[3];
rz(-1.531904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(-1.7588245) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.243467) q[3];
sx q[3];
rz(1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5516978) q[0];
sx q[0];
rz(-2.4015744) q[0];
sx q[0];
rz(-2.9445573) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2824077) q[2];
sx q[2];
rz(-1.7297598) q[2];
sx q[2];
rz(-2.8189916) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41457957) q[1];
sx q[1];
rz(-1.0158744) q[1];
sx q[1];
rz(-1.3347866) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53447978) q[3];
sx q[3];
rz(-1.8090873) q[3];
sx q[3];
rz(-1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.08427944) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(-2.9220667) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052208386) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(2.9242933) q[0];
rz(-0.030933881) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-2.4826179) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4188822) q[0];
sx q[0];
rz(-1.5455149) q[0];
sx q[0];
rz(0.77194571) q[0];
rz(-0.89491567) q[2];
sx q[2];
rz(-0.96940982) q[2];
sx q[2];
rz(-2.2611157) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.79170376) q[1];
sx q[1];
rz(-2.3096497) q[1];
sx q[1];
rz(-1.8658584) q[1];
rz(-pi) q[2];
rz(3.072776) q[3];
sx q[3];
rz(-2.2127082) q[3];
sx q[3];
rz(1.2115692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4387681) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(-2.8213815) q[2];
rz(1.6163588) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76072389) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(2.4587801) q[0];
rz(2.071351) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3037198) q[0];
sx q[0];
rz(-1.2178197) q[0];
sx q[0];
rz(0.2210125) q[0];
x q[1];
rz(2.1754873) q[2];
sx q[2];
rz(-1.7793057) q[2];
sx q[2];
rz(1.0143806) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2883816) q[1];
sx q[1];
rz(-2.2182811) q[1];
sx q[1];
rz(-2.7700469) q[1];
rz(-pi) q[2];
rz(0.50001486) q[3];
sx q[3];
rz(-2.4588636) q[3];
sx q[3];
rz(0.4993712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(2.5749717) q[2];
rz(0.9295272) q[3];
sx q[3];
rz(-2.1598787) q[3];
sx q[3];
rz(-1.7738336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.728445) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(0.79968232) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0737338) q[0];
sx q[0];
rz(-1.6556544) q[0];
sx q[0];
rz(-1.2220864) q[0];
rz(-pi) q[1];
rz(-2.9808447) q[2];
sx q[2];
rz(-1.6511917) q[2];
sx q[2];
rz(0.71842566) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.848222) q[1];
sx q[1];
rz(-1.231041) q[1];
sx q[1];
rz(1.0287702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6639789) q[3];
sx q[3];
rz(-0.16371809) q[3];
sx q[3];
rz(1.9660266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(-2.4463859) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.3643684) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.2790537) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(0.002481133) q[2];
sx q[2];
rz(-0.33591349) q[2];
sx q[2];
rz(-2.9813067) q[2];
rz(-1.4064895) q[3];
sx q[3];
rz(-2.3621029) q[3];
sx q[3];
rz(2.7030871) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
