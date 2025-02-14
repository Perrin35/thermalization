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
rz(-1.7718908) q[0];
sx q[0];
rz(-2.7112024) q[0];
sx q[0];
rz(-2.9769467) q[0];
rz(1.4511664) q[1];
sx q[1];
rz(-2.7632406) q[1];
sx q[1];
rz(-0.70986706) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2458017) q[0];
sx q[0];
rz(-1.4404095) q[0];
sx q[0];
rz(-1.3046645) q[0];
x q[1];
rz(1.1426635) q[2];
sx q[2];
rz(-2.0319418) q[2];
sx q[2];
rz(2.3335339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30395092) q[1];
sx q[1];
rz(-1.6198938) q[1];
sx q[1];
rz(1.9242679) q[1];
rz(-2.6534901) q[3];
sx q[3];
rz(-2.8781366) q[3];
sx q[3];
rz(-1.7149705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5300488) q[2];
sx q[2];
rz(-2.7908466) q[2];
sx q[2];
rz(1.9225527) q[2];
rz(0.48398316) q[3];
sx q[3];
rz(-2.2958906) q[3];
sx q[3];
rz(-1.774196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4265863) q[0];
sx q[0];
rz(-2.8026717) q[0];
sx q[0];
rz(-1.6981079) q[0];
rz(1.8679856) q[1];
sx q[1];
rz(-1.0304291) q[1];
sx q[1];
rz(0.66676203) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50567195) q[0];
sx q[0];
rz(-1.2742325) q[0];
sx q[0];
rz(1.7478329) q[0];
x q[1];
rz(3.0452636) q[2];
sx q[2];
rz(-1.4594635) q[2];
sx q[2];
rz(-1.9736612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8774817) q[1];
sx q[1];
rz(-2.5248233) q[1];
sx q[1];
rz(0.82293226) q[1];
x q[2];
rz(-0.35735766) q[3];
sx q[3];
rz(-0.85225058) q[3];
sx q[3];
rz(2.0524764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.71191177) q[2];
sx q[2];
rz(-1.1744171) q[2];
sx q[2];
rz(-2.5453117) q[2];
rz(-1.0239673) q[3];
sx q[3];
rz(-1.3919316) q[3];
sx q[3];
rz(-1.1202687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.4100818) q[0];
sx q[0];
rz(-0.54786587) q[0];
sx q[0];
rz(-1.4672853) q[0];
rz(3.1027377) q[1];
sx q[1];
rz(-1.7170186) q[1];
sx q[1];
rz(-2.255596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666834) q[0];
sx q[0];
rz(-2.0916478) q[0];
sx q[0];
rz(-1.0731927) q[0];
rz(1.6396693) q[2];
sx q[2];
rz(-2.4703272) q[2];
sx q[2];
rz(1.6948989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13809822) q[1];
sx q[1];
rz(-1.5312863) q[1];
sx q[1];
rz(2.7835151) q[1];
rz(-pi) q[2];
x q[2];
rz(2.233611) q[3];
sx q[3];
rz(-3.040979) q[3];
sx q[3];
rz(-0.20357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2140865) q[2];
sx q[2];
rz(-1.4518041) q[2];
sx q[2];
rz(2.5778594) q[2];
rz(0.8756513) q[3];
sx q[3];
rz(-2.0893658) q[3];
sx q[3];
rz(1.0912857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57986528) q[0];
sx q[0];
rz(-0.848333) q[0];
sx q[0];
rz(-1.0680098) q[0];
rz(-2.7470398) q[1];
sx q[1];
rz(-2.6119472) q[1];
sx q[1];
rz(-1.4716757) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1009586) q[0];
sx q[0];
rz(-1.9143595) q[0];
sx q[0];
rz(-2.6885607) q[0];
rz(-1.2023694) q[2];
sx q[2];
rz(-0.30469182) q[2];
sx q[2];
rz(-2.7161689) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.010943451) q[1];
sx q[1];
rz(-2.2211381) q[1];
sx q[1];
rz(0.54231142) q[1];
rz(0.25217523) q[3];
sx q[3];
rz(-2.7311595) q[3];
sx q[3];
rz(-1.3595734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7431405) q[2];
sx q[2];
rz(-1.8299711) q[2];
sx q[2];
rz(-0.12942806) q[2];
rz(-2.5982924) q[3];
sx q[3];
rz(-2.048251) q[3];
sx q[3];
rz(-1.3885385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1035136) q[0];
sx q[0];
rz(-1.0032126) q[0];
sx q[0];
rz(-0.45573086) q[0];
rz(0.56574667) q[1];
sx q[1];
rz(-2.5028298) q[1];
sx q[1];
rz(-2.9260213) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8889897) q[0];
sx q[0];
rz(-0.71160331) q[0];
sx q[0];
rz(1.3093185) q[0];
rz(-pi) q[1];
rz(2.2224765) q[2];
sx q[2];
rz(-1.0813776) q[2];
sx q[2];
rz(-2.8280743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0255007) q[1];
sx q[1];
rz(-1.3057723) q[1];
sx q[1];
rz(0.071556655) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10978635) q[3];
sx q[3];
rz(-1.854114) q[3];
sx q[3];
rz(2.7840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.7171219) q[2];
sx q[2];
rz(-0.30693808) q[2];
sx q[2];
rz(-1.7090939) q[2];
rz(2.0004382) q[3];
sx q[3];
rz(-1.2204095) q[3];
sx q[3];
rz(-0.46877638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6345217) q[0];
sx q[0];
rz(-0.85352007) q[0];
sx q[0];
rz(-2.5216907) q[0];
rz(1.2076521) q[1];
sx q[1];
rz(-2.4687605) q[1];
sx q[1];
rz(-2.8501453) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47675596) q[0];
sx q[0];
rz(-0.58692044) q[0];
sx q[0];
rz(-1.2715601) q[0];
rz(-pi) q[1];
rz(2.7619128) q[2];
sx q[2];
rz(-2.7951193) q[2];
sx q[2];
rz(2.068067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4174704) q[1];
sx q[1];
rz(-0.63878585) q[1];
sx q[1];
rz(-2.7330509) q[1];
rz(1.5813555) q[3];
sx q[3];
rz(-2.227475) q[3];
sx q[3];
rz(-2.9845723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.42346272) q[2];
sx q[2];
rz(-2.4726157) q[2];
sx q[2];
rz(-2.1675229) q[2];
rz(1.918321) q[3];
sx q[3];
rz(-1.7134824) q[3];
sx q[3];
rz(1.0380925) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85422) q[0];
sx q[0];
rz(-1.9725476) q[0];
sx q[0];
rz(-0.30782345) q[0];
rz(0.27990118) q[1];
sx q[1];
rz(-2.8925536) q[1];
sx q[1];
rz(-1.8797967) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9591589) q[0];
sx q[0];
rz(-1.4366581) q[0];
sx q[0];
rz(0.034237513) q[0];
x q[1];
rz(-1.1818188) q[2];
sx q[2];
rz(-2.2643467) q[2];
sx q[2];
rz(2.3626568) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3615001) q[1];
sx q[1];
rz(-0.29528015) q[1];
sx q[1];
rz(1.5321945) q[1];
rz(-pi) q[2];
rz(2.0359751) q[3];
sx q[3];
rz(-1.5081468) q[3];
sx q[3];
rz(2.986547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6383003) q[2];
sx q[2];
rz(-1.952848) q[2];
sx q[2];
rz(0.80257455) q[2];
rz(-1.5447626) q[3];
sx q[3];
rz(-2.7464726) q[3];
sx q[3];
rz(1.6667268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.06642) q[0];
sx q[0];
rz(-1.5479227) q[0];
sx q[0];
rz(2.5286034) q[0];
rz(-2.0892443) q[1];
sx q[1];
rz(-0.78951013) q[1];
sx q[1];
rz(-1.2552415) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6832812) q[0];
sx q[0];
rz(-1.1363239) q[0];
sx q[0];
rz(2.6660454) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1034715) q[2];
sx q[2];
rz(-1.3824859) q[2];
sx q[2];
rz(1.7323158) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.71129942) q[1];
sx q[1];
rz(-1.0984527) q[1];
sx q[1];
rz(0.76763726) q[1];
rz(2.6735821) q[3];
sx q[3];
rz(-0.29230329) q[3];
sx q[3];
rz(1.4128895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5455948) q[2];
sx q[2];
rz(-0.91384807) q[2];
sx q[2];
rz(2.9929898) q[2];
rz(0.654486) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(1.436208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504836) q[0];
sx q[0];
rz(-1.5330667) q[0];
sx q[0];
rz(-0.0012375687) q[0];
rz(-1.2507863) q[1];
sx q[1];
rz(-0.93235278) q[1];
sx q[1];
rz(0.95157448) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5785268) q[0];
sx q[0];
rz(-2.9349083) q[0];
sx q[0];
rz(-0.61949586) q[0];
x q[1];
rz(-2.6926413) q[2];
sx q[2];
rz(-2.1336974) q[2];
sx q[2];
rz(2.721691) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7504362) q[1];
sx q[1];
rz(-1.5215148) q[1];
sx q[1];
rz(1.5261914) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4755842) q[3];
sx q[3];
rz(-1.1723174) q[3];
sx q[3];
rz(1.0131032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1452267) q[2];
sx q[2];
rz(-1.1684912) q[2];
sx q[2];
rz(1.036198) q[2];
rz(0.869831) q[3];
sx q[3];
rz(-1.6706322) q[3];
sx q[3];
rz(0.74365348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0319801) q[0];
sx q[0];
rz(-0.59015048) q[0];
sx q[0];
rz(-1.1085283) q[0];
rz(2.176586) q[1];
sx q[1];
rz(-1.3771649) q[1];
sx q[1];
rz(-0.56934294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8764502) q[0];
sx q[0];
rz(-2.5776349) q[0];
sx q[0];
rz(-2.7258858) q[0];
x q[1];
rz(-0.37306786) q[2];
sx q[2];
rz(-1.1343435) q[2];
sx q[2];
rz(1.8907569) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.20726897) q[1];
sx q[1];
rz(-2.6187463) q[1];
sx q[1];
rz(-1.5873413) q[1];
rz(-pi) q[2];
rz(-0.17579349) q[3];
sx q[3];
rz(-2.3758278) q[3];
sx q[3];
rz(2.898223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1278594) q[2];
sx q[2];
rz(-2.6467085) q[2];
sx q[2];
rz(0.046517046) q[2];
rz(-2.8401996) q[3];
sx q[3];
rz(-1.2208341) q[3];
sx q[3];
rz(-0.2218328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2409869) q[0];
sx q[0];
rz(-1.2099246) q[0];
sx q[0];
rz(-1.8151617) q[0];
rz(-1.9164512) q[1];
sx q[1];
rz(-0.50347181) q[1];
sx q[1];
rz(1.0540963) q[1];
rz(-1.3780807) q[2];
sx q[2];
rz(-2.1753934) q[2];
sx q[2];
rz(2.8080804) q[2];
rz(-0.71618373) q[3];
sx q[3];
rz(-0.89792218) q[3];
sx q[3];
rz(-0.69178892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
