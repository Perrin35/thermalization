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
rz(-2.5118339) q[0];
sx q[0];
rz(-0.30183733) q[0];
sx q[0];
rz(-1.8344185) q[0];
rz(-0.52840003) q[1];
sx q[1];
rz(-0.83146787) q[1];
sx q[1];
rz(0.45933476) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0889657) q[0];
sx q[0];
rz(-1.3148493) q[0];
sx q[0];
rz(-0.050764485) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1856543) q[2];
sx q[2];
rz(-0.3345851) q[2];
sx q[2];
rz(2.0841887) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3078008) q[1];
sx q[1];
rz(-1.7437309) q[1];
sx q[1];
rz(0.61638919) q[1];
rz(-pi) q[2];
rz(-0.90199043) q[3];
sx q[3];
rz(-1.588527) q[3];
sx q[3];
rz(0.56541857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1186195) q[2];
sx q[2];
rz(-1.3381253) q[2];
sx q[2];
rz(0.55707923) q[2];
rz(0.45595512) q[3];
sx q[3];
rz(-2.3796701) q[3];
sx q[3];
rz(-2.5383811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93422455) q[0];
sx q[0];
rz(-0.71894431) q[0];
sx q[0];
rz(1.7508605) q[0];
rz(-1.3181744) q[1];
sx q[1];
rz(-1.428182) q[1];
sx q[1];
rz(-2.9255829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0681422) q[0];
sx q[0];
rz(-1.2552069) q[0];
sx q[0];
rz(-2.1852949) q[0];
x q[1];
rz(-1.6902709) q[2];
sx q[2];
rz(-1.154197) q[2];
sx q[2];
rz(-2.4753776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60571721) q[1];
sx q[1];
rz(-2.247896) q[1];
sx q[1];
rz(0.89800055) q[1];
rz(-1.029929) q[3];
sx q[3];
rz(-1.3824302) q[3];
sx q[3];
rz(0.68205591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4231437) q[2];
sx q[2];
rz(-2.7687912) q[2];
sx q[2];
rz(-0.049840363) q[2];
rz(2.5981264) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(1.0928833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922888) q[0];
sx q[0];
rz(-1.1019305) q[0];
sx q[0];
rz(-0.9683384) q[0];
rz(-0.020542055) q[1];
sx q[1];
rz(-1.7612709) q[1];
sx q[1];
rz(-1.7113908) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0413301) q[0];
sx q[0];
rz(-1.5484637) q[0];
sx q[0];
rz(-2.5002527) q[0];
x q[1];
rz(1.9400575) q[2];
sx q[2];
rz(-2.103269) q[2];
sx q[2];
rz(-0.42064253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.80950981) q[1];
sx q[1];
rz(-0.44594279) q[1];
sx q[1];
rz(-0.46582241) q[1];
rz(-0.3768206) q[3];
sx q[3];
rz(-1.6018036) q[3];
sx q[3];
rz(-3.0135807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.7057544) q[2];
sx q[2];
rz(-2.4987554) q[2];
sx q[2];
rz(-1.359681) q[2];
rz(-0.5082353) q[3];
sx q[3];
rz(-1.619092) q[3];
sx q[3];
rz(-1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0559167) q[0];
sx q[0];
rz(-0.29656947) q[0];
sx q[0];
rz(-2.1348409) q[0];
rz(2.309917) q[1];
sx q[1];
rz(-2.3995903) q[1];
sx q[1];
rz(1.1078018) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2336034) q[0];
sx q[0];
rz(-2.395326) q[0];
sx q[0];
rz(2.1265592) q[0];
rz(-1.318827) q[2];
sx q[2];
rz(-1.9356079) q[2];
sx q[2];
rz(-3.0598726) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9398168) q[1];
sx q[1];
rz(-0.30758938) q[1];
sx q[1];
rz(-1.5053476) q[1];
rz(-1.8516225) q[3];
sx q[3];
rz(-2.0005595) q[3];
sx q[3];
rz(0.44441477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6882249) q[2];
sx q[2];
rz(-1.0658762) q[2];
sx q[2];
rz(0.31309703) q[2];
rz(0.39737663) q[3];
sx q[3];
rz(-1.4704967) q[3];
sx q[3];
rz(2.9562922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50936407) q[0];
sx q[0];
rz(-0.86357421) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(-0.46982345) q[1];
sx q[1];
rz(-1.3146105) q[1];
sx q[1];
rz(-2.125461) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7454769) q[0];
sx q[0];
rz(-1.1594311) q[0];
sx q[0];
rz(-1.0883254) q[0];
x q[1];
rz(2.8723628) q[2];
sx q[2];
rz(-1.9105663) q[2];
sx q[2];
rz(1.7328078) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.713607) q[1];
sx q[1];
rz(-1.3283037) q[1];
sx q[1];
rz(0.20771435) q[1];
rz(-0.031948745) q[3];
sx q[3];
rz(-1.2476139) q[3];
sx q[3];
rz(-1.972071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4887345) q[2];
sx q[2];
rz(-2.8159339) q[2];
sx q[2];
rz(-0.48750901) q[2];
rz(-0.89753914) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(-1.0645617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2010736) q[0];
sx q[0];
rz(-2.2036393) q[0];
sx q[0];
rz(2.4679389) q[0];
rz(-1.2334088) q[1];
sx q[1];
rz(-2.1234832) q[1];
sx q[1];
rz(2.3755551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2904486) q[0];
sx q[0];
rz(-0.58935114) q[0];
sx q[0];
rz(-1.0714156) q[0];
rz(-pi) q[1];
rz(-0.17579097) q[2];
sx q[2];
rz(-1.5002709) q[2];
sx q[2];
rz(0.27276892) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1704857) q[1];
sx q[1];
rz(-1.0909214) q[1];
sx q[1];
rz(-2.9779153) q[1];
x q[2];
rz(-0.41292878) q[3];
sx q[3];
rz(-2.8097162) q[3];
sx q[3];
rz(1.4982893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.172714) q[2];
sx q[2];
rz(-1.4977027) q[2];
sx q[2];
rz(2.3806351) q[2];
rz(0.40192762) q[3];
sx q[3];
rz(-2.8765078) q[3];
sx q[3];
rz(-0.5886122) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589979) q[0];
sx q[0];
rz(-2.3880279) q[0];
sx q[0];
rz(-1.6931417) q[0];
rz(-2.6407369) q[1];
sx q[1];
rz(-1.8468937) q[1];
sx q[1];
rz(-0.81333152) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7194448) q[0];
sx q[0];
rz(-1.9564183) q[0];
sx q[0];
rz(-1.411012) q[0];
rz(-3.0221239) q[2];
sx q[2];
rz(-1.0257693) q[2];
sx q[2];
rz(-1.8785005) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5593181) q[1];
sx q[1];
rz(-2.6804003) q[1];
sx q[1];
rz(-2.8282437) q[1];
x q[2];
rz(-0.058815033) q[3];
sx q[3];
rz(-1.6928738) q[3];
sx q[3];
rz(-2.5466145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27998754) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(-0.32312265) q[2];
rz(1.4019639) q[3];
sx q[3];
rz(-1.8596545) q[3];
sx q[3];
rz(1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0793656) q[0];
sx q[0];
rz(-1.0492188) q[0];
sx q[0];
rz(2.3754689) q[0];
rz(0.8194204) q[1];
sx q[1];
rz(-1.0304281) q[1];
sx q[1];
rz(1.6105509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19056828) q[0];
sx q[0];
rz(-1.6389209) q[0];
sx q[0];
rz(2.9124252) q[0];
x q[1];
rz(-2.859191) q[2];
sx q[2];
rz(-0.67495433) q[2];
sx q[2];
rz(2.9815428) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8360525) q[1];
sx q[1];
rz(-0.96287268) q[1];
sx q[1];
rz(-2.4120283) q[1];
rz(-pi) q[2];
rz(-2.2977912) q[3];
sx q[3];
rz(-0.42359023) q[3];
sx q[3];
rz(0.18720489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7649585) q[2];
sx q[2];
rz(-2.9464293) q[2];
sx q[2];
rz(-0.39806077) q[2];
rz(2.5352488) q[3];
sx q[3];
rz(-1.5747993) q[3];
sx q[3];
rz(0.91355598) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23676087) q[0];
sx q[0];
rz(-2.8396711) q[0];
sx q[0];
rz(-2.6171369) q[0];
rz(-0.66871387) q[1];
sx q[1];
rz(-1.6122183) q[1];
sx q[1];
rz(-1.761577) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3789744) q[0];
sx q[0];
rz(-1.5272104) q[0];
sx q[0];
rz(0.99751465) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0247714) q[2];
sx q[2];
rz(-0.54440343) q[2];
sx q[2];
rz(-2.4449206) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7126969) q[1];
sx q[1];
rz(-1.98723) q[1];
sx q[1];
rz(2.685999) q[1];
rz(-pi) q[2];
x q[2];
rz(1.947375) q[3];
sx q[3];
rz(-2.0717952) q[3];
sx q[3];
rz(-0.11160103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1606719) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(0.46372947) q[2];
rz(1.7175698) q[3];
sx q[3];
rz(-1.0878891) q[3];
sx q[3];
rz(0.033528479) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5793107) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(0.41807362) q[0];
rz(-0.75830165) q[1];
sx q[1];
rz(-0.98379358) q[1];
sx q[1];
rz(-1.8565348) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4649974) q[0];
sx q[0];
rz(-0.9228068) q[0];
sx q[0];
rz(-1.6654439) q[0];
rz(-pi) q[1];
rz(-1.1056771) q[2];
sx q[2];
rz(-2.4659202) q[2];
sx q[2];
rz(-0.0059520324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.549379) q[1];
sx q[1];
rz(-2.5587808) q[1];
sx q[1];
rz(-2.7874376) q[1];
rz(-pi) q[2];
x q[2];
rz(2.140652) q[3];
sx q[3];
rz(-0.39467282) q[3];
sx q[3];
rz(0.81854023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.16984223) q[2];
sx q[2];
rz(-1.0047793) q[2];
sx q[2];
rz(-0.24610914) q[2];
rz(1.2717815) q[3];
sx q[3];
rz(-1.5747986) q[3];
sx q[3];
rz(-1.7790214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1590189) q[0];
sx q[0];
rz(-1.4763426) q[0];
sx q[0];
rz(2.939298) q[0];
rz(2.5166439) q[1];
sx q[1];
rz(-2.2849871) q[1];
sx q[1];
rz(-2.7402592) q[1];
rz(-0.18904674) q[2];
sx q[2];
rz(-1.3658267) q[2];
sx q[2];
rz(1.7455802) q[2];
rz(2.128849) q[3];
sx q[3];
rz(-2.774917) q[3];
sx q[3];
rz(-2.9493368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
