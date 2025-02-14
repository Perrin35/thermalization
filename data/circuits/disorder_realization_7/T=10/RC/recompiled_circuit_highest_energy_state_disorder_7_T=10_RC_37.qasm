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
rz(1.3697019) q[0];
sx q[0];
rz(-0.43039027) q[0];
sx q[0];
rz(-0.16464591) q[0];
rz(1.4511664) q[1];
sx q[1];
rz(-2.7632406) q[1];
sx q[1];
rz(-0.70986706) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89579099) q[0];
sx q[0];
rz(-1.7011832) q[0];
sx q[0];
rz(-1.3046645) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4998998) q[2];
sx q[2];
rz(-1.1898433) q[2];
sx q[2];
rz(-0.56239201) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8376417) q[1];
sx q[1];
rz(-1.6198938) q[1];
sx q[1];
rz(1.2173247) q[1];
rz(-0.23386896) q[3];
sx q[3];
rz(-1.4483671) q[3];
sx q[3];
rz(-2.8120638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5300488) q[2];
sx q[2];
rz(-0.350746) q[2];
sx q[2];
rz(-1.2190399) q[2];
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
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71500635) q[0];
sx q[0];
rz(-2.8026717) q[0];
sx q[0];
rz(1.6981079) q[0];
rz(-1.8679856) q[1];
sx q[1];
rz(-1.0304291) q[1];
sx q[1];
rz(2.4748306) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979157) q[0];
sx q[0];
rz(-2.7975492) q[0];
sx q[0];
rz(-2.6187569) q[0];
rz(-2.2812541) q[2];
sx q[2];
rz(-2.9945012) q[2];
sx q[2];
rz(1.8837613) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8774817) q[1];
sx q[1];
rz(-0.6167694) q[1];
sx q[1];
rz(2.3186604) q[1];
rz(-pi) q[2];
rz(-0.35735766) q[3];
sx q[3];
rz(-2.2893421) q[3];
sx q[3];
rz(1.0891162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71191177) q[2];
sx q[2];
rz(-1.9671755) q[2];
sx q[2];
rz(2.5453117) q[2];
rz(-2.1176254) q[3];
sx q[3];
rz(-1.3919316) q[3];
sx q[3];
rz(-2.021324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.4100818) q[0];
sx q[0];
rz(-0.54786587) q[0];
sx q[0];
rz(1.6743073) q[0];
rz(0.038854988) q[1];
sx q[1];
rz(-1.424574) q[1];
sx q[1];
rz(0.88599667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9666834) q[0];
sx q[0];
rz(-1.0499448) q[0];
sx q[0];
rz(1.0731927) q[0];
rz(-3.0869834) q[2];
sx q[2];
rz(-2.24018) q[2];
sx q[2];
rz(1.3588248) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3275324) q[1];
sx q[1];
rz(-2.7814354) q[1];
sx q[1];
rz(-0.11231695) q[1];
rz(-pi) q[2];
rz(1.4913849) q[3];
sx q[3];
rz(-1.5089499) q[3];
sx q[3];
rz(-1.1140149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2140865) q[2];
sx q[2];
rz(-1.4518041) q[2];
sx q[2];
rz(-0.56373325) q[2];
rz(-2.2659414) q[3];
sx q[3];
rz(-1.0522269) q[3];
sx q[3];
rz(-1.0912857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57986528) q[0];
sx q[0];
rz(-2.2932597) q[0];
sx q[0];
rz(-2.0735829) q[0];
rz(-0.39455286) q[1];
sx q[1];
rz(-2.6119472) q[1];
sx q[1];
rz(-1.669917) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36763299) q[0];
sx q[0];
rz(-1.1460222) q[0];
sx q[0];
rz(1.9494809) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.285423) q[2];
sx q[2];
rz(-1.6790519) q[2];
sx q[2];
rz(-1.6433604) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.010943451) q[1];
sx q[1];
rz(-0.92045451) q[1];
sx q[1];
rz(-2.5992812) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4626462) q[3];
sx q[3];
rz(-1.9675072) q[3];
sx q[3];
rz(1.0856398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7431405) q[2];
sx q[2];
rz(-1.8299711) q[2];
sx q[2];
rz(0.12942806) q[2];
rz(-0.5433003) q[3];
sx q[3];
rz(-1.0933417) q[3];
sx q[3];
rz(1.7530542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.038079) q[0];
sx q[0];
rz(-2.1383801) q[0];
sx q[0];
rz(-0.45573086) q[0];
rz(2.575846) q[1];
sx q[1];
rz(-0.63876286) q[1];
sx q[1];
rz(0.21557132) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087061398) q[0];
sx q[0];
rz(-2.2534459) q[0];
sx q[0];
rz(-0.21933098) q[0];
x q[1];
rz(0.59026123) q[2];
sx q[2];
rz(-2.1356842) q[2];
sx q[2];
rz(1.5399726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6681155) q[1];
sx q[1];
rz(-1.6398506) q[1];
sx q[1];
rz(-1.8364688) q[1];
x q[2];
rz(1.2858538) q[3];
sx q[3];
rz(-1.4654034) q[3];
sx q[3];
rz(1.8975429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7171219) q[2];
sx q[2];
rz(-2.8346546) q[2];
sx q[2];
rz(1.4324987) q[2];
rz(-2.0004382) q[3];
sx q[3];
rz(-1.9211831) q[3];
sx q[3];
rz(2.6728163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5070709) q[0];
sx q[0];
rz(-0.85352007) q[0];
sx q[0];
rz(-0.61990196) q[0];
rz(1.2076521) q[1];
sx q[1];
rz(-2.4687605) q[1];
sx q[1];
rz(-2.8501453) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3454764) q[0];
sx q[0];
rz(-1.7347851) q[0];
sx q[0];
rz(-1.0046474) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4377777) q[2];
sx q[2];
rz(-1.2499193) q[2];
sx q[2];
rz(1.6668298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7241223) q[1];
sx q[1];
rz(-0.63878585) q[1];
sx q[1];
rz(-0.40854172) q[1];
x q[2];
rz(0.013697946) q[3];
sx q[3];
rz(-2.4848416) q[3];
sx q[3];
rz(2.9672772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42346272) q[2];
sx q[2];
rz(-2.4726157) q[2];
sx q[2];
rz(-2.1675229) q[2];
rz(-1.918321) q[3];
sx q[3];
rz(-1.4281102) q[3];
sx q[3];
rz(1.0380925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85422) q[0];
sx q[0];
rz(-1.9725476) q[0];
sx q[0];
rz(-2.8337692) q[0];
rz(-2.8616915) q[1];
sx q[1];
rz(-2.8925536) q[1];
sx q[1];
rz(-1.8797967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7578106) q[0];
sx q[0];
rz(-1.5368665) q[0];
sx q[0];
rz(-1.4365804) q[0];
rz(2.4096411) q[2];
sx q[2];
rz(-1.2748655) q[2];
sx q[2];
rz(2.0934806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9692291) q[1];
sx q[1];
rz(-1.5595655) q[1];
sx q[1];
rz(-1.8658691) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1056176) q[3];
sx q[3];
rz(-1.5081468) q[3];
sx q[3];
rz(0.1550457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5032924) q[2];
sx q[2];
rz(-1.952848) q[2];
sx q[2];
rz(-0.80257455) q[2];
rz(-1.59683) q[3];
sx q[3];
rz(-2.7464726) q[3];
sx q[3];
rz(-1.6667268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751727) q[0];
sx q[0];
rz(-1.5936699) q[0];
sx q[0];
rz(2.5286034) q[0];
rz(2.0892443) q[1];
sx q[1];
rz(-2.3520825) q[1];
sx q[1];
rz(-1.2552415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57277623) q[0];
sx q[0];
rz(-2.5089087) q[0];
sx q[0];
rz(-2.3494472) q[0];
rz(0.21771149) q[2];
sx q[2];
rz(-1.048511) q[2];
sx q[2];
rz(3.0899855) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.44574983) q[1];
sx q[1];
rz(-0.90409213) q[1];
sx q[1];
rz(0.95335828) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26237285) q[3];
sx q[3];
rz(-1.4404357) q[3];
sx q[3];
rz(-2.8487132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5959979) q[2];
sx q[2];
rz(-2.2277446) q[2];
sx q[2];
rz(-2.9929898) q[2];
rz(0.654486) q[3];
sx q[3];
rz(-1.196967) q[3];
sx q[3];
rz(1.436208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-0.29110903) q[0];
sx q[0];
rz(-1.608526) q[0];
sx q[0];
rz(-0.0012375687) q[0];
rz(-1.2507863) q[1];
sx q[1];
rz(-2.2092399) q[1];
sx q[1];
rz(-0.95157448) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5785268) q[0];
sx q[0];
rz(-2.9349083) q[0];
sx q[0];
rz(0.61949586) q[0];
rz(-pi) q[1];
rz(-2.6926413) q[2];
sx q[2];
rz(-1.0078953) q[2];
sx q[2];
rz(-2.721691) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7504362) q[1];
sx q[1];
rz(-1.6200778) q[1];
sx q[1];
rz(1.5261914) q[1];
rz(-2.9195027) q[3];
sx q[3];
rz(-0.40910334) q[3];
sx q[3];
rz(2.3698185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1452267) q[2];
sx q[2];
rz(-1.1684912) q[2];
sx q[2];
rz(2.1053947) q[2];
rz(2.2717617) q[3];
sx q[3];
rz(-1.6706322) q[3];
sx q[3];
rz(2.3979392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1096126) q[0];
sx q[0];
rz(-0.59015048) q[0];
sx q[0];
rz(2.0330644) q[0];
rz(0.96500665) q[1];
sx q[1];
rz(-1.7644278) q[1];
sx q[1];
rz(-0.56934294) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78380221) q[0];
sx q[0];
rz(-1.0598425) q[0];
sx q[0];
rz(-1.3207256) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7685248) q[2];
sx q[2];
rz(-2.0072492) q[2];
sx q[2];
rz(1.8907569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.20726897) q[1];
sx q[1];
rz(-0.52284635) q[1];
sx q[1];
rz(-1.5542514) q[1];
x q[2];
rz(-2.9657992) q[3];
sx q[3];
rz(-2.3758278) q[3];
sx q[3];
rz(-2.898223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0137332) q[2];
sx q[2];
rz(-2.6467085) q[2];
sx q[2];
rz(-3.0950756) q[2];
rz(2.8401996) q[3];
sx q[3];
rz(-1.9207585) q[3];
sx q[3];
rz(-0.2218328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90060577) q[0];
sx q[0];
rz(-1.931668) q[0];
sx q[0];
rz(1.3264309) q[0];
rz(-1.9164512) q[1];
sx q[1];
rz(-0.50347181) q[1];
sx q[1];
rz(1.0540963) q[1];
rz(2.8711748) q[2];
sx q[2];
rz(-2.5106988) q[2];
sx q[2];
rz(-0.0028263447) q[2];
rz(-0.75792652) q[3];
sx q[3];
rz(-1.0316385) q[3];
sx q[3];
rz(1.3760174) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
