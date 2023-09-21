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
rz(-2.8074582) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(-1.2184719) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96555644) q[0];
sx q[0];
rz(-2.3210225) q[0];
sx q[0];
rz(-1.5216989) q[0];
x q[1];
rz(0.46484868) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(2.6236617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7114746) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(0.75573604) q[1];
rz(-pi) q[2];
rz(0.86745947) q[3];
sx q[3];
rz(-0.28775035) q[3];
sx q[3];
rz(0.61144637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-0.66453385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.739025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70799202) q[0];
sx q[0];
rz(-3.1118244) q[0];
sx q[0];
rz(0.13694163) q[0];
rz(2.5820929) q[2];
sx q[2];
rz(-1.7569949) q[2];
sx q[2];
rz(-0.49636832) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84239292) q[1];
sx q[1];
rz(-2.0626915) q[1];
sx q[1];
rz(1.1115587) q[1];
rz(-pi) q[2];
rz(-1.4757231) q[3];
sx q[3];
rz(-0.94082309) q[3];
sx q[3];
rz(-0.41439393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.3690255) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(0.40107045) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6221878) q[0];
sx q[0];
rz(-1.6739968) q[0];
sx q[0];
rz(1.2445356) q[0];
rz(-0.6443278) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(0.49567859) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.229278) q[1];
sx q[1];
rz(-2.4114128) q[1];
sx q[1];
rz(3.1320523) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0173666) q[3];
sx q[3];
rz(-0.43101573) q[3];
sx q[3];
rz(-0.22180804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(-0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(-2.8881883) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69617535) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(-2.1131383) q[0];
rz(-pi) q[1];
x q[1];
rz(0.01732145) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(-2.3522365) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.464444) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(1.4518751) q[1];
x q[2];
rz(-3.0074189) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(-1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040745) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(0.0054365693) q[0];
x q[1];
rz(-0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(-0.82816154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3921567) q[1];
sx q[1];
rz(-1.4538308) q[1];
sx q[1];
rz(-2.3163296) q[1];
x q[2];
rz(3.0466988) q[3];
sx q[3];
rz(-2.0695544) q[3];
sx q[3];
rz(2.9961078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(2.9120973) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5620419) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29612449) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(2.092093) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74395545) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(-0.26847408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8333203) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(-2.1182548) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(-3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(1.1118719) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3355013) q[0];
sx q[0];
rz(-1.4772692) q[0];
sx q[0];
rz(-2.0557687) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9260336) q[2];
sx q[2];
rz(-2.07395) q[2];
sx q[2];
rz(-2.9687198) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8060311) q[1];
sx q[1];
rz(-2.8415488) q[1];
sx q[1];
rz(0.26144822) q[1];
x q[2];
rz(-2.0734378) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(2.2964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-2.4107966) q[0];
rz(0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98159957) q[0];
sx q[0];
rz(-0.85708517) q[0];
sx q[0];
rz(0.1937565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63379143) q[2];
sx q[2];
rz(-2.3049424) q[2];
sx q[2];
rz(-1.447669) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6983812) q[1];
sx q[1];
rz(-1.4581919) q[1];
sx q[1];
rz(2.5618782) q[1];
rz(-pi) q[2];
rz(2.6492277) q[3];
sx q[3];
rz(-2.1015321) q[3];
sx q[3];
rz(2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-0.63684741) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.586097) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(0.019502217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6396128) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(-1.7777068) q[0];
x q[1];
rz(-1.3703913) q[2];
sx q[2];
rz(-1.7192171) q[2];
sx q[2];
rz(-0.60147775) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2724185) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(-1.8723349) q[1];
x q[2];
rz(-1.8832302) q[3];
sx q[3];
rz(-1.7282681) q[3];
sx q[3];
rz(-1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(-1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432817) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(1.9932995) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023037993) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(1.5492357) q[0];
rz(-pi) q[1];
rz(-0.63126385) q[2];
sx q[2];
rz(-1.751465) q[2];
sx q[2];
rz(2.0618912) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1710098) q[1];
sx q[1];
rz(-0.90921558) q[1];
sx q[1];
rz(-0.2172825) q[1];
rz(-2.9989472) q[3];
sx q[3];
rz(-0.35257617) q[3];
sx q[3];
rz(-0.25776097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0317595) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(-2.6504436) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
rz(3.1046042) q[3];
sx q[3];
rz(-2.130641) q[3];
sx q[3];
rz(-0.11161042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
