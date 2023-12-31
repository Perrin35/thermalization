OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1887814) q[0];
sx q[0];
rz(-0.45438284) q[0];
sx q[0];
rz(2.7812468) q[0];
rz(-0.25720815) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(0.53127015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60018051) q[1];
sx q[1];
rz(-1.0804847) q[1];
sx q[1];
rz(0.4835101) q[1];
rz(2.1336018) q[3];
sx q[3];
rz(-1.6562914) q[3];
sx q[3];
rz(2.8443955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(2.2757754) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(3.1153733) q[0];
rz(-1.5401309) q[1];
sx q[1];
rz(-1.5988348) q[1];
sx q[1];
rz(0.96347934) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.547077) q[0];
sx q[0];
rz(-1.5684959) q[0];
sx q[0];
rz(0.95675795) q[0];
rz(-pi) q[1];
rz(1.0019148) q[2];
sx q[2];
rz(-1.2895673) q[2];
sx q[2];
rz(2.1525454) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0986833) q[1];
sx q[1];
rz(-2.5334362) q[1];
sx q[1];
rz(-2.152918) q[1];
x q[2];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.8027455) q[3];
sx q[3];
rz(1.9333145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-3.0070686) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(2.1988595) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117675) q[0];
sx q[0];
rz(-0.38914248) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(0.55999666) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6277498) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(1.9613128) q[0];
rz(-1.5921002) q[2];
sx q[2];
rz(-1.2676123) q[2];
sx q[2];
rz(1.6456749) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3126038) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(2.7752084) q[1];
rz(-pi) q[2];
rz(0.47645724) q[3];
sx q[3];
rz(-1.1343079) q[3];
sx q[3];
rz(-0.95526327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75227633) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(-2.9690572) q[2];
rz(-2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(1.0579695) q[3];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1383706) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.8428615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772229) q[0];
sx q[0];
rz(-1.7958613) q[0];
sx q[0];
rz(-1.0610915) q[0];
x q[1];
rz(-0.0012871731) q[2];
sx q[2];
rz(-2.3372071) q[2];
sx q[2];
rz(-0.019891642) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5609834) q[1];
sx q[1];
rz(-2.1316075) q[1];
sx q[1];
rz(0.22052712) q[1];
rz(-pi) q[2];
rz(1.2508568) q[3];
sx q[3];
rz(-2.7971929) q[3];
sx q[3];
rz(-0.26564769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6056885) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.6872905) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6032747) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(2.8447661) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083858) q[0];
sx q[0];
rz(-2.1413681) q[0];
sx q[0];
rz(-1.5242566) q[0];
rz(-1.7237687) q[2];
sx q[2];
rz(-2.3077871) q[2];
sx q[2];
rz(2.134915) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0592812) q[1];
sx q[1];
rz(-1.1059522) q[1];
sx q[1];
rz(2.6962198) q[1];
x q[2];
rz(-3.0389901) q[3];
sx q[3];
rz(-0.71628621) q[3];
sx q[3];
rz(-0.41075452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.8313235) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(-1.1522419) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(0.31744441) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-2.3902067) q[0];
rz(-1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-0.60633916) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060121814) q[0];
sx q[0];
rz(-1.524964) q[0];
sx q[0];
rz(-1.5541374) q[0];
x q[1];
rz(-0.41646429) q[2];
sx q[2];
rz(-0.93566862) q[2];
sx q[2];
rz(3.0481899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50748435) q[1];
sx q[1];
rz(-2.1180696) q[1];
sx q[1];
rz(-0.21168153) q[1];
rz(-pi) q[2];
rz(-1.353225) q[3];
sx q[3];
rz(-1.6440344) q[3];
sx q[3];
rz(-0.53461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(-0.50810057) q[0];
rz(1.5628901) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9433141) q[0];
sx q[0];
rz(-2.4662848) q[0];
sx q[0];
rz(2.6576463) q[0];
rz(-pi) q[1];
rz(1.5825282) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(-2.8649883) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8040647) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(-0.40562628) q[1];
rz(-pi) q[2];
rz(0.64738691) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(-0.90741531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.7283758) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(-3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(1.3887127) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213288) q[0];
sx q[0];
rz(-1.5409894) q[0];
sx q[0];
rz(-2.1304312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65865626) q[2];
sx q[2];
rz(-0.63630644) q[2];
sx q[2];
rz(-3.0898526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82909225) q[1];
sx q[1];
rz(-1.478985) q[1];
sx q[1];
rz(-1.5593668) q[1];
rz(2.1771168) q[3];
sx q[3];
rz(-1.8165605) q[3];
sx q[3];
rz(2.5442459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.7154988) q[0];
rz(-0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-0.55823278) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8919864) q[0];
sx q[0];
rz(-2.5476646) q[0];
sx q[0];
rz(2.3114165) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6782645) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(1.4449643) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3840752) q[1];
sx q[1];
rz(-1.5069403) q[1];
sx q[1];
rz(-2.3149895) q[1];
x q[2];
rz(0.34096034) q[3];
sx q[3];
rz(-2.6302359) q[3];
sx q[3];
rz(-0.90358678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(-0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(0.82540712) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(0.5232946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494203) q[0];
sx q[0];
rz(-1.6329375) q[0];
sx q[0];
rz(0.6092351) q[0];
x q[1];
rz(0.9988437) q[2];
sx q[2];
rz(-1.6076644) q[2];
sx q[2];
rz(-1.8429304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1112422) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(0.72279795) q[1];
rz(-pi) q[2];
rz(-2.8966122) q[3];
sx q[3];
rz(-1.3462726) q[3];
sx q[3];
rz(2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.837073) q[2];
sx q[2];
rz(-2.4261116) q[2];
sx q[2];
rz(-0.26930299) q[2];
rz(2.6473911) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.3080296) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(0.79694637) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
