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
rz(1.310362) q[1];
sx q[1];
rz(5.789776) q[1];
sx q[1];
rz(8.9763666) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7767186) q[0];
sx q[0];
rz(-1.5176511) q[0];
sx q[0];
rz(2.4928635) q[0];
rz(-pi) q[1];
rz(-2.7712819) q[2];
sx q[2];
rz(-2.8680621) q[2];
sx q[2];
rz(0.47765484) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55576217) q[1];
sx q[1];
rz(-2.1735648) q[1];
sx q[1];
rz(0.65048154) q[1];
rz(-3.0541522) q[3];
sx q[3];
rz(-2.6387847) q[3];
sx q[3];
rz(1.8016165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51556921) q[2];
sx q[2];
rz(-1.3884576) q[2];
sx q[2];
rz(0.62421978) q[2];
rz(-0.37912399) q[3];
sx q[3];
rz(-2.1513394) q[3];
sx q[3];
rz(-0.24851255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666331) q[0];
sx q[0];
rz(-0.81460726) q[0];
sx q[0];
rz(-2.1061184) q[0];
rz(1.2738719) q[1];
sx q[1];
rz(-0.73718166) q[1];
sx q[1];
rz(2.6587291) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55601701) q[0];
sx q[0];
rz(-3.034798) q[0];
sx q[0];
rz(0.44321816) q[0];
rz(-pi) q[1];
rz(0.16629433) q[2];
sx q[2];
rz(-1.6876843) q[2];
sx q[2];
rz(2.6720195) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.00512) q[1];
sx q[1];
rz(-2.0025952) q[1];
sx q[1];
rz(3.1312607) q[1];
x q[2];
rz(-2.8206536) q[3];
sx q[3];
rz(-1.648099) q[3];
sx q[3];
rz(-2.2805205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.024293385) q[2];
sx q[2];
rz(-1.6397986) q[2];
sx q[2];
rz(-1.8710322) q[2];
rz(-2.471586) q[3];
sx q[3];
rz(-0.62205258) q[3];
sx q[3];
rz(-0.89699927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73827493) q[0];
sx q[0];
rz(-0.44602317) q[0];
sx q[0];
rz(0.73905149) q[0];
rz(-2.1801379) q[1];
sx q[1];
rz(-2.8196204) q[1];
sx q[1];
rz(-0.75698537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3630545) q[0];
sx q[0];
rz(-0.86150384) q[0];
sx q[0];
rz(0.34687931) q[0];
x q[1];
rz(-0.69125533) q[2];
sx q[2];
rz(-1.8658085) q[2];
sx q[2];
rz(1.6694836) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9582102) q[1];
sx q[1];
rz(-2.498495) q[1];
sx q[1];
rz(2.9166168) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0717172) q[3];
sx q[3];
rz(-1.2638076) q[3];
sx q[3];
rz(-1.2941293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6387393) q[2];
sx q[2];
rz(-2.2245378) q[2];
sx q[2];
rz(-2.4364831) q[2];
rz(2.8201568) q[3];
sx q[3];
rz(-2.1240081) q[3];
sx q[3];
rz(-0.41079918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0099156378) q[0];
sx q[0];
rz(-2.3035045) q[0];
sx q[0];
rz(-1.2257082) q[0];
rz(1.907584) q[1];
sx q[1];
rz(-2.1261413) q[1];
sx q[1];
rz(3.0753678) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28540885) q[0];
sx q[0];
rz(-1.3254306) q[0];
sx q[0];
rz(-3.0491327) q[0];
x q[1];
rz(-0.80133665) q[2];
sx q[2];
rz(-1.7883915) q[2];
sx q[2];
rz(2.269553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.855763) q[1];
sx q[1];
rz(-2.0501839) q[1];
sx q[1];
rz(-2.781032) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31853557) q[3];
sx q[3];
rz(-2.7580166) q[3];
sx q[3];
rz(0.84795241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0951198) q[2];
sx q[2];
rz(-2.1472223) q[2];
sx q[2];
rz(-2.7009916) q[2];
rz(1.0062086) q[3];
sx q[3];
rz(-1.3738084) q[3];
sx q[3];
rz(0.83427507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42136583) q[0];
sx q[0];
rz(-0.81291968) q[0];
sx q[0];
rz(-1.6356069) q[0];
rz(1.5274564) q[1];
sx q[1];
rz(-1.2200049) q[1];
sx q[1];
rz(-1.45586) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86690534) q[0];
sx q[0];
rz(-2.2969249) q[0];
sx q[0];
rz(-0.19498904) q[0];
rz(1.6873932) q[2];
sx q[2];
rz(-1.3449838) q[2];
sx q[2];
rz(-1.1979529) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33259049) q[1];
sx q[1];
rz(-2.2223516) q[1];
sx q[1];
rz(1.5452191) q[1];
rz(-pi) q[2];
rz(2.9221228) q[3];
sx q[3];
rz(-2.5220036) q[3];
sx q[3];
rz(2.0989024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.27318925) q[2];
sx q[2];
rz(-1.578178) q[2];
sx q[2];
rz(0.91147649) q[2];
rz(2.2677926) q[3];
sx q[3];
rz(-2.3995212) q[3];
sx q[3];
rz(-3.1349283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8580496) q[0];
sx q[0];
rz(-0.25074211) q[0];
sx q[0];
rz(0.13033303) q[0];
rz(0.48769543) q[1];
sx q[1];
rz(-2.6002488) q[1];
sx q[1];
rz(0.74388751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34602964) q[0];
sx q[0];
rz(-1.5209805) q[0];
sx q[0];
rz(1.6532142) q[0];
rz(3.0024372) q[2];
sx q[2];
rz(-0.94991131) q[2];
sx q[2];
rz(0.54229743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2762017) q[1];
sx q[1];
rz(-1.8542093) q[1];
sx q[1];
rz(2.1332425) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33244953) q[3];
sx q[3];
rz(-2.6327193) q[3];
sx q[3];
rz(-0.32989855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0222212) q[2];
sx q[2];
rz(-2.2890942) q[2];
sx q[2];
rz(0.96088299) q[2];
rz(-0.39572257) q[3];
sx q[3];
rz(-1.8053677) q[3];
sx q[3];
rz(0.1161639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6845067) q[0];
sx q[0];
rz(-2.0261903) q[0];
sx q[0];
rz(-1.047026) q[0];
rz(-2.537435) q[1];
sx q[1];
rz(-1.8641169) q[1];
sx q[1];
rz(0.39271694) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0020802845) q[0];
sx q[0];
rz(-0.89186984) q[0];
sx q[0];
rz(2.9154791) q[0];
rz(-pi) q[1];
rz(-0.057633295) q[2];
sx q[2];
rz(-0.95256348) q[2];
sx q[2];
rz(1.0317486) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3369357) q[1];
sx q[1];
rz(-1.5976904) q[1];
sx q[1];
rz(0.27534816) q[1];
rz(-pi) q[2];
rz(-2.4684577) q[3];
sx q[3];
rz(-2.7302448) q[3];
sx q[3];
rz(-2.213221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0543694) q[2];
sx q[2];
rz(-1.1944218) q[2];
sx q[2];
rz(1.8325904) q[2];
rz(-1.7878112) q[3];
sx q[3];
rz(-1.196922) q[3];
sx q[3];
rz(0.27031171) q[3];
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
rz(0.29276174) q[0];
sx q[0];
rz(-1.4563541) q[0];
sx q[0];
rz(-2.8344179) q[0];
rz(1.3245026) q[1];
sx q[1];
rz(-2.7565286) q[1];
sx q[1];
rz(-2.2408392) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4221638) q[0];
sx q[0];
rz(-1.3086645) q[0];
sx q[0];
rz(-1.5882701) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45489725) q[2];
sx q[2];
rz(-2.1730246) q[2];
sx q[2];
rz(-1.0554316) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72462725) q[1];
sx q[1];
rz(-0.46675439) q[1];
sx q[1];
rz(-2.3326567) q[1];
rz(-2.9489904) q[3];
sx q[3];
rz(-0.64036548) q[3];
sx q[3];
rz(0.34512025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33425346) q[2];
sx q[2];
rz(-0.05143493) q[2];
sx q[2];
rz(2.5267498) q[2];
rz(-2.0452512) q[3];
sx q[3];
rz(-2.5494826) q[3];
sx q[3];
rz(-0.69826564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39356247) q[0];
sx q[0];
rz(-2.5108971) q[0];
sx q[0];
rz(0.50931859) q[0];
rz(2.9604984) q[1];
sx q[1];
rz(-1.4053922) q[1];
sx q[1];
rz(-0.96955713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80705331) q[0];
sx q[0];
rz(-1.4167804) q[0];
sx q[0];
rz(-2.3898983) q[0];
rz(3.0907057) q[2];
sx q[2];
rz(-2.0774088) q[2];
sx q[2];
rz(-1.2801054) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0694094) q[1];
sx q[1];
rz(-1.9095712) q[1];
sx q[1];
rz(1.4563926) q[1];
rz(-pi) q[2];
rz(1.9593616) q[3];
sx q[3];
rz(-2.3747184) q[3];
sx q[3];
rz(-2.9076613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2481044) q[2];
sx q[2];
rz(-0.23427811) q[2];
sx q[2];
rz(1.0815557) q[2];
rz(2.9099416) q[3];
sx q[3];
rz(-1.3737498) q[3];
sx q[3];
rz(-0.20802465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3807826) q[0];
sx q[0];
rz(-2.91687) q[0];
sx q[0];
rz(-2.494452) q[0];
rz(-0.45516792) q[1];
sx q[1];
rz(-1.4249233) q[1];
sx q[1];
rz(-2.3796577) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1200346) q[0];
sx q[0];
rz(-1.2823449) q[0];
sx q[0];
rz(0.7201654) q[0];
x q[1];
rz(2.2438887) q[2];
sx q[2];
rz(-0.87397777) q[2];
sx q[2];
rz(-0.76329939) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0207723) q[1];
sx q[1];
rz(-1.5365531) q[1];
sx q[1];
rz(-2.4832151) q[1];
x q[2];
rz(-1.1033789) q[3];
sx q[3];
rz(-1.5776792) q[3];
sx q[3];
rz(-0.21871834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0882988) q[2];
sx q[2];
rz(-2.9247354) q[2];
sx q[2];
rz(2.2376412) q[2];
rz(1.5229185) q[3];
sx q[3];
rz(-2.121033) q[3];
sx q[3];
rz(1.1618377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.224613) q[0];
sx q[0];
rz(-1.9539178) q[0];
sx q[0];
rz(1.7896347) q[0];
rz(3.0147973) q[1];
sx q[1];
rz(-0.8225816) q[1];
sx q[1];
rz(-2.2093028) q[1];
rz(-2.8895072) q[2];
sx q[2];
rz(-2.0093976) q[2];
sx q[2];
rz(-1.2653399) q[2];
rz(-0.0062777304) q[3];
sx q[3];
rz(-2.7115887) q[3];
sx q[3];
rz(-0.030621519) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
