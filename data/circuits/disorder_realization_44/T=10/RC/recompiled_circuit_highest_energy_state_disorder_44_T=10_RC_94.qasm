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
rz(-2.0353844) q[0];
sx q[0];
rz(-0.68618965) q[0];
sx q[0];
rz(-1.1535147) q[0];
rz(-1.8312307) q[1];
sx q[1];
rz(-2.6481833) q[1];
sx q[1];
rz(-2.6931813) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8657239) q[0];
sx q[0];
rz(-2.4910036) q[0];
sx q[0];
rz(-0.087819789) q[0];
x q[1];
rz(-2.885778) q[2];
sx q[2];
rz(-1.4728776) q[2];
sx q[2];
rz(-2.4061587) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.55576217) q[1];
sx q[1];
rz(-2.1735648) q[1];
sx q[1];
rz(2.4911111) q[1];
rz(3.0541522) q[3];
sx q[3];
rz(-0.50280792) q[3];
sx q[3];
rz(-1.3399762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(-0.62421978) q[2];
rz(2.7624687) q[3];
sx q[3];
rz(-2.1513394) q[3];
sx q[3];
rz(-0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17495951) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(-2.1061184) q[0];
rz(-1.2738719) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(0.48286352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5855756) q[0];
sx q[0];
rz(-3.034798) q[0];
sx q[0];
rz(2.6983745) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16629433) q[2];
sx q[2];
rz(-1.4539084) q[2];
sx q[2];
rz(0.4695732) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9804364) q[1];
sx q[1];
rz(-2.709678) q[1];
sx q[1];
rz(1.5483787) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9008174) q[3];
sx q[3];
rz(-0.32980439) q[3];
sx q[3];
rz(-2.2036116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1172993) q[2];
sx q[2];
rz(-1.5017941) q[2];
sx q[2];
rz(-1.2705605) q[2];
rz(-0.67000669) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(-2.2445934) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(-2.4025412) q[0];
rz(2.1801379) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(-2.3846073) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3630545) q[0];
sx q[0];
rz(-0.86150384) q[0];
sx q[0];
rz(-0.34687931) q[0];
rz(-0.44481014) q[2];
sx q[2];
rz(-2.399657) q[2];
sx q[2];
rz(-2.7049989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18338246) q[1];
sx q[1];
rz(-2.498495) q[1];
sx q[1];
rz(-0.22497589) q[1];
rz(-pi) q[2];
rz(1.3540165) q[3];
sx q[3];
rz(-0.31459537) q[3];
sx q[3];
rz(1.6198688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5028533) q[2];
sx q[2];
rz(-0.91705489) q[2];
sx q[2];
rz(-2.4364831) q[2];
rz(2.8201568) q[3];
sx q[3];
rz(-1.0175846) q[3];
sx q[3];
rz(-2.7307935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(3.131677) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(1.2257082) q[0];
rz(1.2340087) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(-3.0753678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8787251) q[0];
sx q[0];
rz(-1.6604794) q[0];
sx q[0];
rz(-1.32442) q[0];
rz(0.80133665) q[2];
sx q[2];
rz(-1.7883915) q[2];
sx q[2];
rz(0.8720397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5410616) q[1];
sx q[1];
rz(-2.5503073) q[1];
sx q[1];
rz(-0.97452428) q[1];
rz(-2.7755967) q[3];
sx q[3];
rz(-1.4533236) q[3];
sx q[3];
rz(-2.7155196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0464728) q[2];
sx q[2];
rz(-0.9943704) q[2];
sx q[2];
rz(-2.7009916) q[2];
rz(2.1353841) q[3];
sx q[3];
rz(-1.7677842) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42136583) q[0];
sx q[0];
rz(-2.328673) q[0];
sx q[0];
rz(1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.2200049) q[1];
sx q[1];
rz(-1.45586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5680926) q[0];
sx q[0];
rz(-1.4254036) q[0];
sx q[0];
rz(2.306434) q[0];
rz(-pi) q[1];
rz(-0.46874502) q[2];
sx q[2];
rz(-0.25368099) q[2];
sx q[2];
rz(1.4616592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2537186) q[1];
sx q[1];
rz(-1.5504596) q[1];
sx q[1];
rz(2.4898797) q[1];
rz(-2.5334356) q[3];
sx q[3];
rz(-1.4440315) q[3];
sx q[3];
rz(-2.7931282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27318925) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(-2.2301162) q[2];
rz(-0.8738001) q[3];
sx q[3];
rz(-0.74207145) q[3];
sx q[3];
rz(-0.0066643683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28354302) q[0];
sx q[0];
rz(-2.8908505) q[0];
sx q[0];
rz(-0.13033303) q[0];
rz(-0.48769543) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(-0.74388751) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9209393) q[0];
sx q[0];
rz(-1.488481) q[0];
sx q[0];
rz(-0.049985188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3792453) q[2];
sx q[2];
rz(-2.5073176) q[2];
sx q[2];
rz(2.3630362) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6106657) q[1];
sx q[1];
rz(-2.1082889) q[1];
sx q[1];
rz(0.33157886) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7508932) q[3];
sx q[3];
rz(-1.0922179) q[3];
sx q[3];
rz(-0.046600051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0222212) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(-2.7458701) q[3];
sx q[3];
rz(-1.3362249) q[3];
sx q[3];
rz(0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.457086) q[0];
sx q[0];
rz(-2.0261903) q[0];
sx q[0];
rz(2.0945666) q[0];
rz(-0.60415769) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(-0.39271694) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0020802845) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(-2.9154791) q[0];
x q[1];
rz(-0.95177841) q[2];
sx q[2];
rz(-1.6177532) q[2];
sx q[2];
rz(-2.6359735) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8046569) q[1];
sx q[1];
rz(-1.5439022) q[1];
sx q[1];
rz(-2.8662445) q[1];
x q[2];
rz(0.32870558) q[3];
sx q[3];
rz(-1.3188601) q[3];
sx q[3];
rz(-3.1303034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0543694) q[2];
sx q[2];
rz(-1.9471709) q[2];
sx q[2];
rz(1.8325904) q[2];
rz(1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(-0.27031171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8488309) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(2.8344179) q[0];
rz(-1.81709) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-2.2408392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3548292) q[0];
sx q[0];
rz(-2.8788924) q[0];
sx q[0];
rz(-3.0765615) q[0];
rz(-pi) q[1];
rz(-1.0020761) q[2];
sx q[2];
rz(-0.7373215) q[2];
sx q[2];
rz(-1.3739283) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5985721) q[1];
sx q[1];
rz(-1.9024339) q[1];
sx q[1];
rz(2.806862) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9489904) q[3];
sx q[3];
rz(-0.64036548) q[3];
sx q[3];
rz(-2.7964724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33425346) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(-0.61484289) q[2];
rz(-2.0452512) q[3];
sx q[3];
rz(-0.59211007) q[3];
sx q[3];
rz(0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480302) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(-2.6322741) q[0];
rz(0.18109426) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(0.96955713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2354202) q[0];
sx q[0];
rz(-2.3114822) q[0];
sx q[0];
rz(1.3614015) q[0];
x q[1];
rz(-1.0636341) q[2];
sx q[2];
rz(-1.5263057) q[2];
sx q[2];
rz(-0.26598334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.5395551) q[1];
sx q[1];
rz(-1.4629211) q[1];
sx q[1];
rz(-2.8007568) q[1];
rz(-pi) q[2];
rz(0.35004444) q[3];
sx q[3];
rz(-0.87348233) q[3];
sx q[3];
rz(-2.3907135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2481044) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(-1.0815557) q[2];
rz(2.9099416) q[3];
sx q[3];
rz(-1.7678429) q[3];
sx q[3];
rz(-2.933568) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3807826) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(2.494452) q[0];
rz(2.6864247) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(-2.3796577) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3477701) q[0];
sx q[0];
rz(-0.88632727) q[0];
sx q[0];
rz(-1.1948299) q[0];
rz(2.2438887) q[2];
sx q[2];
rz(-2.2676149) q[2];
sx q[2];
rz(0.76329939) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.40574351) q[1];
sx q[1];
rz(-2.4824576) q[1];
sx q[1];
rz(3.0856611) q[1];
x q[2];
rz(-3.1338828) q[3];
sx q[3];
rz(-2.0382018) q[3];
sx q[3];
rz(-1.3555519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0882988) q[2];
sx q[2];
rz(-0.21685728) q[2];
sx q[2];
rz(-2.2376412) q[2];
rz(1.6186742) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(1.9797549) q[3];
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
rz(-0.91697964) q[0];
sx q[0];
rz(-1.1876748) q[0];
sx q[0];
rz(-1.351958) q[0];
rz(0.12679535) q[1];
sx q[1];
rz(-2.3190111) q[1];
sx q[1];
rz(0.93228985) q[1];
rz(-0.2520855) q[2];
sx q[2];
rz(-1.132195) q[2];
sx q[2];
rz(1.8762527) q[2];
rz(1.5736754) q[3];
sx q[3];
rz(-2.0007912) q[3];
sx q[3];
rz(3.1040647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
