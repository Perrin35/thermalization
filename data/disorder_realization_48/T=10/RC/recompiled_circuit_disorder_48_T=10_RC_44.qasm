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
rz(-1.2184719) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2479808) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(-0.052608842) q[0];
x q[1];
rz(2.676744) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(0.51793098) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.76518607) q[1];
sx q[1];
rz(-2.0284488) q[1];
sx q[1];
rz(-0.54995579) q[1];
x q[2];
rz(-2.2741332) q[3];
sx q[3];
rz(-2.8538423) q[3];
sx q[3];
rz(2.5301463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(2.4770588) q[3];
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
rz(-2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.739025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57099045) q[0];
sx q[0];
rz(-1.6002858) q[0];
sx q[0];
rz(-1.5667314) q[0];
rz(-pi) q[1];
rz(1.7895133) q[2];
sx q[2];
rz(-2.1195076) q[2];
sx q[2];
rz(0.95900853) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1837511) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(2.6028231) q[1];
rz(-2.5094633) q[3];
sx q[3];
rz(-1.4940133) q[3];
sx q[3];
rz(2.0413105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7188321) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(-1.4136219) q[0];
rz(2.6638022) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.385752) q[0];
sx q[0];
rz(-0.34163654) q[0];
sx q[0];
rz(-1.2582448) q[0];
rz(-pi) q[1];
rz(0.24984078) q[2];
sx q[2];
rz(-2.4820538) q[2];
sx q[2];
rz(1.8674873) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9123147) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(0.0095403949) q[1];
rz(-pi) q[2];
rz(1.1242261) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(-2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(-0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.697631) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-2.4855721) q[1];
sx q[1];
rz(2.8881883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69617535) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(-1.0284543) q[0];
rz(-2.0565815) q[2];
sx q[2];
rz(-1.5554785) q[2];
sx q[2];
rz(0.77335301) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6654012) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(0.16237662) q[1];
rz(-pi) q[2];
rz(3.0074189) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(-1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(0.45809349) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3375181) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-3.1361561) q[0];
rz(0.96180054) q[2];
sx q[2];
rz(-0.70448175) q[2];
sx q[2];
rz(1.5722164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0702857) q[1];
sx q[1];
rz(-2.3100393) q[1];
sx q[1];
rz(0.15858312) q[1];
rz(-1.3985653) q[3];
sx q[3];
rz(-2.6346364) q[3];
sx q[3];
rz(-0.34190049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.3389448) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(-0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.2671635) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29612449) q[0];
sx q[0];
rz(-0.86588973) q[0];
sx q[0];
rz(-1.0494997) q[0];
x q[1];
rz(0.74395545) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(-2.8731186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3082723) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(2.1182548) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24354981) q[3];
sx q[3];
rz(-2.0381513) q[3];
sx q[3];
rz(-1.5598284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(-1.0725347) q[0];
rz(-0.5468927) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(2.0297208) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5892964) q[0];
sx q[0];
rz(-0.49320212) q[0];
sx q[0];
rz(1.769355) q[0];
rz(-pi) q[1];
x q[1];
rz(1.215559) q[2];
sx q[2];
rz(-1.0676427) q[2];
sx q[2];
rz(0.1728729) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.062473) q[1];
sx q[1];
rz(-1.8603431) q[1];
sx q[1];
rz(1.6505961) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0681549) q[3];
sx q[3];
rz(-1.916269) q[3];
sx q[3];
rz(-0.84514602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(1.973935) q[2];
rz(1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(1.4060098) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.3370445) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98159957) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(0.1937565) q[0];
x q[1];
rz(-2.1515498) q[2];
sx q[2];
rz(-2.2120737) q[2];
sx q[2];
rz(2.5255447) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6983812) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(-0.57971445) q[1];
rz(-pi) q[2];
rz(-0.98324361) q[3];
sx q[3];
rz(-1.9907111) q[3];
sx q[3];
rz(-2.0411232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-0.019502217) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6396128) q[0];
sx q[0];
rz(-1.6507971) q[0];
sx q[0];
rz(-1.7777068) q[0];
x q[1];
rz(1.3703913) q[2];
sx q[2];
rz(-1.7192171) q[2];
sx q[2];
rz(-2.5401149) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86917415) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(1.2692578) q[1];
x q[2];
rz(2.0476258) q[3];
sx q[3];
rz(-2.7928824) q[3];
sx q[3];
rz(2.204012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(-1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(1.1482931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1295706) q[0];
sx q[0];
rz(-0.14229933) q[0];
sx q[0];
rz(2.9905031) q[0];
x q[1];
rz(2.5103288) q[2];
sx q[2];
rz(-1.3901276) q[2];
sx q[2];
rz(1.0797015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53459586) q[1];
sx q[1];
rz(-1.3998704) q[1];
sx q[1];
rz(-2.2439438) q[1];
x q[2];
rz(-0.14264543) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(-0.25776097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2293573) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(-1.9671494) q[2];
sx q[2];
rz(-2.195993) q[2];
sx q[2];
rz(3.1396951) q[2];
rz(-1.6297324) q[3];
sx q[3];
rz(-0.56093506) q[3];
sx q[3];
rz(-0.18118071) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
