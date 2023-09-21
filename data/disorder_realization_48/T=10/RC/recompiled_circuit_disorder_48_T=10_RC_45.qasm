OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(0.33413449) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1760362) q[0];
sx q[0];
rz(-0.82057014) q[0];
sx q[0];
rz(1.6198938) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6115509) q[2];
sx q[2];
rz(-1.1062804) q[2];
sx q[2];
rz(-2.0704616) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3764066) q[1];
sx q[1];
rz(-1.1131439) q[1];
sx q[1];
rz(-2.5916369) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3487885) q[3];
sx q[3];
rz(-1.3862002) q[3];
sx q[3];
rz(0.27664646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1419066) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(-0.02948972) q[0];
rz(-pi) q[1];
rz(-1.7895133) q[2];
sx q[2];
rz(-2.1195076) q[2];
sx q[2];
rz(-0.95900853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1837511) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(-0.5387696) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0121147) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(-0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(-0.31769162) q[2];
rz(-2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(2.3247705) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0553592) q[0];
sx q[0];
rz(-1.8952574) q[0];
sx q[0];
rz(-0.1089036) q[0];
rz(-0.24984078) q[2];
sx q[2];
rz(-2.4820538) q[2];
sx q[2];
rz(1.2741054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.89951) q[1];
sx q[1];
rz(-0.84065719) q[1];
sx q[1];
rz(1.5793369) q[1];
x q[2];
rz(-1.1242261) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(2.3613789) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-0.25340432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4454173) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(2.1131383) q[0];
x q[1];
rz(1.6035945) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4761915) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(-0.16237662) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0074189) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(0.45809349) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8100909) q[0];
sx q[0];
rz(-2.0134263) q[0];
sx q[0];
rz(-1.5733733) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1796218) q[2];
sx q[2];
rz(-1.95032) q[2];
sx q[2];
rz(0.48987197) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.304368) q[1];
sx q[1];
rz(-2.3886884) q[1];
sx q[1];
rz(1.7423082) q[1];
rz(-pi) q[2];
rz(-1.3985653) q[3];
sx q[3];
rz(-2.6346364) q[3];
sx q[3];
rz(-0.34190049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(-0.22949533) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5620419) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.2671635) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6309109) q[0];
sx q[0];
rz(-1.9598538) q[0];
sx q[0];
rz(-2.3657777) q[0];
rz(-pi) q[1];
rz(2.9847449) q[2];
sx q[2];
rz(-0.75011293) q[2];
sx q[2];
rz(-1.4175121) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3082723) q[1];
sx q[1];
rz(-0.38494021) q[1];
sx q[1];
rz(2.1182548) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24354981) q[3];
sx q[3];
rz(-1.1034414) q[3];
sx q[3];
rz(-1.5817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(1.1118719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3355013) q[0];
sx q[0];
rz(-1.6643235) q[0];
sx q[0];
rz(-2.0557687) q[0];
rz(-pi) q[1];
rz(-2.6107437) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(1.9206778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.48549451) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(-0.2904201) q[1];
rz(-pi) q[2];
rz(2.2124347) q[3];
sx q[3];
rz(-2.540179) q[3];
sx q[3];
rz(1.8638368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.1676577) q[2];
rz(1.6052823) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-0.73079601) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98159957) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(0.1937565) q[0];
x q[1];
rz(-2.412699) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(-0.58065562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29754408) q[1];
sx q[1];
rz(-2.5522759) q[1];
sx q[1];
rz(2.938016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(-2.0411232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(-2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(3.1220904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0855904) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(0.081736728) q[0];
rz(1.3703913) q[2];
sx q[2];
rz(-1.4223756) q[2];
sx q[2];
rz(-0.60147775) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86917415) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(-1.2692578) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0939668) q[3];
sx q[3];
rz(-0.34871021) q[3];
sx q[3];
rz(-2.204012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(-2.7539339) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(0.48450255) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(1.1482931) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185547) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(1.592357) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7933153) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(-2.7811108) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3155568) q[1];
sx q[1];
rz(-2.450374) q[1];
sx q[1];
rz(-1.8408937) q[1];
rz(-pi) q[2];
rz(-0.14264543) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(2.8838317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(0.036988463) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];