OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463319) q[0];
sx q[0];
rz(-1.594461) q[0];
sx q[0];
rz(-0.46121009) q[0];
rz(-pi) q[1];
rz(1.2167591) q[2];
sx q[2];
rz(-2.3166222) q[2];
sx q[2];
rz(1.1302469) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2156148) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(1.2222626) q[1];
rz(-pi) q[2];
rz(0.10856467) q[3];
sx q[3];
rz(-1.8130455) q[3];
sx q[3];
rz(3.0755175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5550845) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(-2.5449975) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.2260431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2873043) q[0];
sx q[0];
rz(-2.0137631) q[0];
sx q[0];
rz(1.7914377) q[0];
x q[1];
rz(2.3184899) q[2];
sx q[2];
rz(-2.1777993) q[2];
sx q[2];
rz(-2.5603103) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.023149816) q[1];
sx q[1];
rz(-2.9217302) q[1];
sx q[1];
rz(1.5727732) q[1];
rz(1.5561043) q[3];
sx q[3];
rz(-0.5726632) q[3];
sx q[3];
rz(0.0088012561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.9689955) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(-0.80667574) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(-1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817292) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(-1.249041) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-2.1121315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91988504) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(-2.1973781) q[0];
rz(2.4918409) q[2];
sx q[2];
rz(-1.8012815) q[2];
sx q[2];
rz(2.8284555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0975115) q[1];
sx q[1];
rz(-1.71002) q[1];
sx q[1];
rz(0.15686762) q[1];
rz(-pi) q[2];
rz(0.73996468) q[3];
sx q[3];
rz(-1.660534) q[3];
sx q[3];
rz(-1.9268074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(-2.6339445) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(-1.1631789) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5959651) q[0];
rz(1.0428628) q[1];
sx q[1];
rz(-1.1735801) q[1];
sx q[1];
rz(1.625659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3665873) q[0];
sx q[0];
rz(-0.89582743) q[0];
sx q[0];
rz(1.7540027) q[0];
rz(1.4342986) q[2];
sx q[2];
rz(-0.56328661) q[2];
sx q[2];
rz(2.0493281) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.792946) q[1];
sx q[1];
rz(-0.9325087) q[1];
sx q[1];
rz(1.9116886) q[1];
rz(-2.8193072) q[3];
sx q[3];
rz(-0.61338131) q[3];
sx q[3];
rz(-1.8005288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.5198583) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-2.246726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8653523) q[0];
sx q[0];
rz(-1.1603174) q[0];
sx q[0];
rz(2.4549237) q[0];
rz(-pi) q[1];
rz(2.1331482) q[2];
sx q[2];
rz(-0.92698669) q[2];
sx q[2];
rz(1.2959727) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6579073) q[1];
sx q[1];
rz(-0.92805082) q[1];
sx q[1];
rz(0.55773463) q[1];
x q[2];
rz(1.5785061) q[3];
sx q[3];
rz(-1.2843411) q[3];
sx q[3];
rz(0.25659284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52508369) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(1.9990702) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.3775795) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87100055) q[0];
sx q[0];
rz(-0.56068476) q[0];
sx q[0];
rz(2.2339348) q[0];
rz(1.6443406) q[2];
sx q[2];
rz(-1.2300756) q[2];
sx q[2];
rz(-1.3982915) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28593674) q[1];
sx q[1];
rz(-0.93614139) q[1];
sx q[1];
rz(-3.082285) q[1];
rz(2.898071) q[3];
sx q[3];
rz(-0.72391073) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12748195) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.4200462) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(0.15596095) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6316846) q[0];
sx q[0];
rz(-1.317306) q[0];
sx q[0];
rz(1.3939199) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1440802) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(-1.13525) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5676463) q[1];
sx q[1];
rz(-1.624755) q[1];
sx q[1];
rz(-1.7458564) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74219269) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.069313958) q[2];
sx q[2];
rz(-0.57702714) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1919365) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(-2.8930194) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0524806) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(0.46778932) q[0];
rz(-pi) q[1];
rz(-0.82069223) q[2];
sx q[2];
rz(-2.4889915) q[2];
sx q[2];
rz(-1.5936268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8057115) q[1];
sx q[1];
rz(-0.98166775) q[1];
sx q[1];
rz(-2.5475403) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.076898889) q[3];
sx q[3];
rz(-1.7401164) q[3];
sx q[3];
rz(2.3016735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8481855) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(-3.0548813) q[2];
rz(-2.6596206) q[3];
sx q[3];
rz(-2.635699) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(0.28276643) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.3185906) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5134207) q[0];
sx q[0];
rz(-1.6512198) q[0];
sx q[0];
rz(-1.1286939) q[0];
rz(-0.65638541) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(0.93751794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2576037) q[1];
sx q[1];
rz(-0.89234967) q[1];
sx q[1];
rz(-0.18737327) q[1];
rz(0.41140822) q[3];
sx q[3];
rz(-2.3477926) q[3];
sx q[3];
rz(-2.0421162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(-2.7424157) q[2];
rz(-0.88360751) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76319641) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(-1.0653161) q[1];
sx q[1];
rz(-2.1677446) q[1];
sx q[1];
rz(-1.8803966) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0364089) q[0];
sx q[0];
rz(-1.9793708) q[0];
sx q[0];
rz(-1.4559792) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59826675) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(-2.604904) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1120158) q[1];
sx q[1];
rz(-0.50260168) q[1];
sx q[1];
rz(-3.0081248) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68393771) q[3];
sx q[3];
rz(-1.6947019) q[3];
sx q[3];
rz(0.77139664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(-1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8284843) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(-2.5333511) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-0.026253168) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(2.5504997) q[3];
sx q[3];
rz(-1.4553087) q[3];
sx q[3];
rz(1.6607264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
