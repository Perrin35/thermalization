OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(1.9231208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63872913) q[0];
sx q[0];
rz(-1.5348866) q[0];
sx q[0];
rz(-0.7508276) q[0];
rz(-pi) q[1];
rz(1.6115509) q[2];
sx q[2];
rz(-2.0353122) q[2];
sx q[2];
rz(-1.0711311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54108809) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(-2.0946676) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86745947) q[3];
sx q[3];
rz(-2.8538423) q[3];
sx q[3];
rz(-2.5301463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(2.4770588) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.4025677) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419066) q[0];
sx q[0];
rz(-1.5748595) q[0];
sx q[0];
rz(0.02948972) q[0];
rz(-pi) q[1];
rz(1.7895133) q[2];
sx q[2];
rz(-1.0220851) q[2];
sx q[2];
rz(-0.95900853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2991997) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(1.1115587) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0121147) q[3];
sx q[3];
rz(-0.63614142) q[3];
sx q[3];
rz(0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-0.31769162) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51940489) q[0];
sx q[0];
rz(-1.4675958) q[0];
sx q[0];
rz(-1.897057) q[0];
rz(2.4972649) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(-0.49567859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.89951) q[1];
sx q[1];
rz(-0.84065719) q[1];
sx q[1];
rz(-1.5793369) q[1];
rz(0.19604711) q[3];
sx q[3];
rz(-1.9571597) q[3];
sx q[3];
rz(2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5473189) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(0.55580124) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-0.25340432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70595104) q[0];
sx q[0];
rz(-2.5905529) q[0];
sx q[0];
rz(-1.7680697) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5379982) q[2];
sx q[2];
rz(-0.48600733) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6771486) q[1];
sx q[1];
rz(-2.1975937) q[1];
sx q[1];
rz(1.6897175) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4686618) q[3];
sx q[3];
rz(-2.2176952) q[3];
sx q[3];
rz(1.1598066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-0.78732642) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-0.45809349) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23819085) q[0];
sx q[0];
rz(-1.5684677) q[0];
sx q[0];
rz(-0.44263126) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96180054) q[2];
sx q[2];
rz(-2.4371109) q[2];
sx q[2];
rz(-1.5722164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.304368) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(-1.7423082) q[1];
x q[2];
rz(3.0466988) q[3];
sx q[3];
rz(-2.0695544) q[3];
sx q[3];
rz(2.9961078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(0.0028006639) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(-1.3389448) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(-1.8744291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7125268) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(-2.6119786) q[0];
x q[1];
rz(-2.3976372) q[2];
sx q[2];
rz(-1.4641054) q[2];
sx q[2];
rz(-0.26847408) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2515182) q[1];
sx q[1];
rz(-1.8971844) q[1];
sx q[1];
rz(-2.9337487) q[1];
x q[2];
rz(-2.0503067) q[3];
sx q[3];
rz(-1.7877842) q[3];
sx q[3];
rz(-3.019141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(-1.8976338) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5892964) q[0];
sx q[0];
rz(-0.49320212) q[0];
sx q[0];
rz(1.3722377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53084897) q[2];
sx q[2];
rz(-1.8804272) q[2];
sx q[2];
rz(1.2209148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8060311) q[1];
sx q[1];
rz(-0.30004382) q[1];
sx q[1];
rz(-0.26144822) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0681549) q[3];
sx q[3];
rz(-1.2253237) q[3];
sx q[3];
rz(-2.2964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(0.73079601) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(0.75497595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98159957) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(2.9478361) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99004284) q[2];
sx q[2];
rz(-2.2120737) q[2];
sx q[2];
rz(2.5255447) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8440486) q[1];
sx q[1];
rz(-0.58931671) q[1];
sx q[1];
rz(2.938016) q[1];
rz(-0.98324361) q[3];
sx q[3];
rz(-1.9907111) q[3];
sx q[3];
rz(1.1004694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(0.63684741) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(1.138858) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-3.1220904) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50197983) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(1.7777068) q[0];
rz(0.15140622) q[2];
sx q[2];
rz(-1.7689686) q[2];
sx q[2];
rz(-0.99934794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.49669493) q[1];
sx q[1];
rz(-0.95458191) q[1];
sx q[1];
rz(-2.9195021) q[1];
x q[2];
rz(-0.1653413) q[3];
sx q[3];
rz(-1.2623566) q[3];
sx q[3];
rz(-2.7066018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(1.1482931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023037993) q[0];
sx q[0];
rz(-1.7114637) q[0];
sx q[0];
rz(-1.5492357) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8414367) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(-0.73210994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82603589) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(1.300699) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(-1.9625361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(0.45583367) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(2.6504436) q[2];
sx q[2];
rz(-0.72577234) q[2];
sx q[2];
rz(2.5189248) q[2];
rz(1.6297324) q[3];
sx q[3];
rz(-2.5806576) q[3];
sx q[3];
rz(2.9604119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];