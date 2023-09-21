OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2258423) q[0];
sx q[0];
rz(-0.031210829) q[0];
sx q[0];
rz(2.6565235) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(5.8689868) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.30487) q[0];
sx q[0];
rz(-0.49855907) q[0];
sx q[0];
rz(2.3303633) q[0];
rz(-1.4235052) q[2];
sx q[2];
rz(-2.5669332) q[2];
sx q[2];
rz(1.6342083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1297076) q[1];
sx q[1];
rz(-1.1457448) q[1];
sx q[1];
rz(0.070406291) q[1];
x q[2];
rz(-1.4685417) q[3];
sx q[3];
rz(-1.3073321) q[3];
sx q[3];
rz(1.2276358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9238613) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(1.2565553) q[3];
sx q[3];
rz(-0.50285181) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(0.1698499) q[0];
rz(-2.4376712) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(2.6020715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9233421) q[0];
sx q[0];
rz(-1.5241511) q[0];
sx q[0];
rz(1.1211066) q[0];
rz(-pi) q[1];
rz(0.45692921) q[2];
sx q[2];
rz(-2.3655342) q[2];
sx q[2];
rz(-1.6844144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19903781) q[1];
sx q[1];
rz(-0.53831646) q[1];
sx q[1];
rz(1.1953137) q[1];
rz(-pi) q[2];
rz(0.86124729) q[3];
sx q[3];
rz(-1.1303139) q[3];
sx q[3];
rz(1.7064106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.66723055) q[2];
sx q[2];
rz(-1.903406) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(-1.2438141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77984017) q[0];
sx q[0];
rz(-3.0267974) q[0];
sx q[0];
rz(-0.4483805) q[0];
rz(1.7547296) q[1];
sx q[1];
rz(-1.153839) q[1];
sx q[1];
rz(-2.8853436) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.019779531) q[0];
sx q[0];
rz(-2.461713) q[0];
sx q[0];
rz(2.7133184) q[0];
rz(-pi) q[1];
rz(0.82661144) q[2];
sx q[2];
rz(-0.6140784) q[2];
sx q[2];
rz(2.5342864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5970522) q[1];
sx q[1];
rz(-1.8982732) q[1];
sx q[1];
rz(0.87045963) q[1];
x q[2];
rz(-2.933421) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(2.4908623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8024575) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(1.908196) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(-0.25948778) q[0];
rz(1.9909987) q[1];
sx q[1];
rz(-1.78803) q[1];
sx q[1];
rz(-2.4096699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5607802) q[0];
sx q[0];
rz(-0.73045759) q[0];
sx q[0];
rz(-2.3985582) q[0];
rz(0.133693) q[2];
sx q[2];
rz(-0.72172726) q[2];
sx q[2];
rz(-1.114811) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96958292) q[1];
sx q[1];
rz(-1.3130377) q[1];
sx q[1];
rz(1.7715363) q[1];
x q[2];
rz(-1.6963523) q[3];
sx q[3];
rz(-1.351965) q[3];
sx q[3];
rz(0.20727508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4449473) q[2];
sx q[2];
rz(-1.2735294) q[2];
sx q[2];
rz(2.990492) q[2];
rz(-0.54667306) q[3];
sx q[3];
rz(-2.0918545) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(1.8792101) q[0];
rz(1.4683912) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(-0.79777065) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3482462) q[0];
sx q[0];
rz(-3.0214546) q[0];
sx q[0];
rz(-1.5915464) q[0];
rz(-pi) q[1];
x q[1];
rz(1.351379) q[2];
sx q[2];
rz(-1.2741538) q[2];
sx q[2];
rz(-0.23362939) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.080575374) q[1];
sx q[1];
rz(-1.4772381) q[1];
sx q[1];
rz(-2.1993125) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2410774) q[3];
sx q[3];
rz(-1.1493249) q[3];
sx q[3];
rz(2.4333405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(2.8395555) q[2];
rz(1.9942412) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7323332) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.9858032) q[0];
rz(-2.0571158) q[1];
sx q[1];
rz(-0.98025727) q[1];
sx q[1];
rz(3.0715122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018774059) q[0];
sx q[0];
rz(-0.2304603) q[0];
sx q[0];
rz(-0.83968681) q[0];
rz(-pi) q[1];
rz(0.55307936) q[2];
sx q[2];
rz(-2.2773909) q[2];
sx q[2];
rz(1.7318219) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49680432) q[1];
sx q[1];
rz(-1.4092688) q[1];
sx q[1];
rz(0.5715538) q[1];
x q[2];
rz(1.4962247) q[3];
sx q[3];
rz(-1.6061022) q[3];
sx q[3];
rz(1.1535742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(-0.62136674) q[2];
rz(-1.7012043) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(2.281718) q[0];
rz(-1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(-0.0079356114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36695776) q[0];
sx q[0];
rz(-1.3152221) q[0];
sx q[0];
rz(2.5711683) q[0];
x q[1];
rz(-0.87848778) q[2];
sx q[2];
rz(-1.2307067) q[2];
sx q[2];
rz(-1.0483339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6565894) q[1];
sx q[1];
rz(-2.2409229) q[1];
sx q[1];
rz(2.2336002) q[1];
x q[2];
rz(-2.8273724) q[3];
sx q[3];
rz(-1.946297) q[3];
sx q[3];
rz(-2.7995031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4456711) q[2];
sx q[2];
rz(-1.2233223) q[2];
sx q[2];
rz(-0.41637862) q[2];
rz(1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-2.2369475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1798582) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(-0.36488786) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84425981) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(0.25015932) q[0];
rz(-pi) q[1];
rz(2.084923) q[2];
sx q[2];
rz(-0.80386111) q[2];
sx q[2];
rz(0.67509292) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.9099721) q[1];
sx q[1];
rz(-0.93712229) q[1];
sx q[1];
rz(-0.85068591) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3650465) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(-2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6442948) q[2];
sx q[2];
rz(-0.50540322) q[2];
sx q[2];
rz(-2.0765182) q[2];
rz(0.30125695) q[3];
sx q[3];
rz(-1.732429) q[3];
sx q[3];
rz(1.5208972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(2.705943) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-2.6711617) q[1];
sx q[1];
rz(-0.41697821) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57396736) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(-0.91659878) q[0];
rz(-pi) q[1];
x q[1];
rz(2.142799) q[2];
sx q[2];
rz(-2.1527094) q[2];
sx q[2];
rz(-1.7172161) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0028487) q[1];
sx q[1];
rz(-1.7094304) q[1];
sx q[1];
rz(2.0744051) q[1];
rz(-pi) q[2];
rz(-3.0005089) q[3];
sx q[3];
rz(-1.699563) q[3];
sx q[3];
rz(1.5878549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0372662) q[2];
sx q[2];
rz(-1.6486847) q[2];
sx q[2];
rz(2.9157675) q[2];
rz(0.2078235) q[3];
sx q[3];
rz(-2.4184629) q[3];
sx q[3];
rz(2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.726783) q[0];
sx q[0];
rz(-2.8746958) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(0.95364755) q[1];
sx q[1];
rz(-1.8667659) q[1];
sx q[1];
rz(1.3226002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37227092) q[0];
sx q[0];
rz(-1.1936545) q[0];
sx q[0];
rz(-3.0110714) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1216713) q[2];
sx q[2];
rz(-1.3438864) q[2];
sx q[2];
rz(-0.34044468) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0304071) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(2.2305957) q[1];
rz(2.9880452) q[3];
sx q[3];
rz(-0.921075) q[3];
sx q[3];
rz(-0.49680199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3502729) q[2];
sx q[2];
rz(-1.046317) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(-1.919205) q[3];
sx q[3];
rz(-2.9639769) q[3];
sx q[3];
rz(0.55541742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9983457) q[0];
sx q[0];
rz(-2.2194942) q[0];
sx q[0];
rz(1.0794328) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-1.580711) q[2];
sx q[2];
rz(-2.0799939) q[2];
sx q[2];
rz(-1.6891198) q[2];
rz(-1.8018467) q[3];
sx q[3];
rz(-1.4828724) q[3];
sx q[3];
rz(2.8433269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];