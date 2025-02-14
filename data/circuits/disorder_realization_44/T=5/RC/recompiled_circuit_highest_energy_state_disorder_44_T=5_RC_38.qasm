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
rz(2.9889838) q[0];
sx q[0];
rz(-1.9419365) q[0];
sx q[0];
rz(-2.8359523) q[0];
rz(-2.8513554) q[1];
sx q[1];
rz(-1.4479535) q[1];
sx q[1];
rz(2.1374968) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5619156) q[0];
sx q[0];
rz(-1.5439057) q[0];
sx q[0];
rz(0.057232522) q[0];
rz(-2.2260336) q[2];
sx q[2];
rz(-2.8311074) q[2];
sx q[2];
rz(0.99570292) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7028382) q[1];
sx q[1];
rz(-1.251529) q[1];
sx q[1];
rz(-1.1478018) q[1];
x q[2];
rz(1.6131931) q[3];
sx q[3];
rz(-1.0323845) q[3];
sx q[3];
rz(-1.1841033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4940138) q[2];
sx q[2];
rz(-0.61654377) q[2];
sx q[2];
rz(0.54440633) q[2];
rz(2.7529649) q[3];
sx q[3];
rz(-1.9184687) q[3];
sx q[3];
rz(2.2012034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87043864) q[0];
sx q[0];
rz(-2.1156023) q[0];
sx q[0];
rz(2.0197268) q[0];
rz(1.059996) q[1];
sx q[1];
rz(-0.50299877) q[1];
sx q[1];
rz(-0.88567919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87410418) q[0];
sx q[0];
rz(-0.42952202) q[0];
sx q[0];
rz(-2.3699479) q[0];
x q[1];
rz(-0.8203542) q[2];
sx q[2];
rz(-2.1848896) q[2];
sx q[2];
rz(2.6439553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52351511) q[1];
sx q[1];
rz(-1.6517868) q[1];
sx q[1];
rz(0.063132719) q[1];
rz(-pi) q[2];
rz(0.84195413) q[3];
sx q[3];
rz(-2.9837787) q[3];
sx q[3];
rz(2.006881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9649967) q[2];
sx q[2];
rz(-1.4623888) q[2];
sx q[2];
rz(-0.029335955) q[2];
rz(0.36597478) q[3];
sx q[3];
rz(-0.47593203) q[3];
sx q[3];
rz(2.9717145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40816864) q[0];
sx q[0];
rz(-1.3329196) q[0];
sx q[0];
rz(-0.15342203) q[0];
rz(-2.9452005) q[1];
sx q[1];
rz(-1.6729665) q[1];
sx q[1];
rz(0.82082716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0401413) q[0];
sx q[0];
rz(-0.87803221) q[0];
sx q[0];
rz(0.67130868) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48604934) q[2];
sx q[2];
rz(-2.1450248) q[2];
sx q[2];
rz(-0.42178139) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7559636) q[1];
sx q[1];
rz(-1.5308994) q[1];
sx q[1];
rz(2.7261158) q[1];
rz(-pi) q[2];
rz(-0.36661212) q[3];
sx q[3];
rz(-1.5248581) q[3];
sx q[3];
rz(2.1076155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11744943) q[2];
sx q[2];
rz(-1.3261745) q[2];
sx q[2];
rz(2.2085025) q[2];
rz(-0.34705958) q[3];
sx q[3];
rz(-1.1091899) q[3];
sx q[3];
rz(0.6412653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1953122) q[0];
sx q[0];
rz(-1.1271789) q[0];
sx q[0];
rz(-2.0830925) q[0];
rz(-0.095254101) q[1];
sx q[1];
rz(-1.9922099) q[1];
sx q[1];
rz(2.9375295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.090074) q[0];
sx q[0];
rz(-1.7885889) q[0];
sx q[0];
rz(2.9831332) q[0];
x q[1];
rz(0.80486678) q[2];
sx q[2];
rz(-0.61986065) q[2];
sx q[2];
rz(-1.0484753) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82641027) q[1];
sx q[1];
rz(-2.1433804) q[1];
sx q[1];
rz(-1.3316283) q[1];
rz(-2.4434632) q[3];
sx q[3];
rz(-1.4120201) q[3];
sx q[3];
rz(0.2793437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.05693398) q[2];
sx q[2];
rz(-1.6997507) q[2];
sx q[2];
rz(1.320768) q[2];
rz(-2.1386347) q[3];
sx q[3];
rz(-1.2489677) q[3];
sx q[3];
rz(-0.5654208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4249307) q[0];
sx q[0];
rz(-1.2053763) q[0];
sx q[0];
rz(-0.64111125) q[0];
rz(2.5504316) q[1];
sx q[1];
rz(-1.7008275) q[1];
sx q[1];
rz(1.444918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.648429) q[0];
sx q[0];
rz(-0.38721353) q[0];
sx q[0];
rz(-3.0680543) q[0];
rz(0.21331738) q[2];
sx q[2];
rz(-2.2498395) q[2];
sx q[2];
rz(-0.17903331) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7406297) q[1];
sx q[1];
rz(-1.1006121) q[1];
sx q[1];
rz(-1.2931254) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8376546) q[3];
sx q[3];
rz(-1.0660488) q[3];
sx q[3];
rz(-0.42298906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.558202) q[2];
sx q[2];
rz(-2.148197) q[2];
sx q[2];
rz(0.85303419) q[2];
rz(2.9663626) q[3];
sx q[3];
rz(-1.3220738) q[3];
sx q[3];
rz(0.62293735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.919642) q[0];
sx q[0];
rz(-2.3355244) q[0];
sx q[0];
rz(-2.6958418) q[0];
rz(2.1469965) q[1];
sx q[1];
rz(-2.1371195) q[1];
sx q[1];
rz(1.4882784) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5588682) q[0];
sx q[0];
rz(-1.7877401) q[0];
sx q[0];
rz(0.72175755) q[0];
rz(1.6449261) q[2];
sx q[2];
rz(-1.526381) q[2];
sx q[2];
rz(-1.4786468) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6116029) q[1];
sx q[1];
rz(-1.7771136) q[1];
sx q[1];
rz(-1.8782658) q[1];
x q[2];
rz(-0.078468389) q[3];
sx q[3];
rz(-1.8967617) q[3];
sx q[3];
rz(-2.7808444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.970924) q[2];
sx q[2];
rz(-1.2430151) q[2];
sx q[2];
rz(-1.0260014) q[2];
rz(-2.3435727) q[3];
sx q[3];
rz(-2.3012216) q[3];
sx q[3];
rz(-2.87288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2949424) q[0];
sx q[0];
rz(-1.3947399) q[0];
sx q[0];
rz(-0.41644874) q[0];
rz(1.396748) q[1];
sx q[1];
rz(-1.6114707) q[1];
sx q[1];
rz(1.4769311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0234982) q[0];
sx q[0];
rz(-2.5865285) q[0];
sx q[0];
rz(1.5007625) q[0];
rz(-pi) q[1];
rz(2.8139122) q[2];
sx q[2];
rz(-1.9679356) q[2];
sx q[2];
rz(-0.59573345) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1394964) q[1];
sx q[1];
rz(-1.1285892) q[1];
sx q[1];
rz(0.79331974) q[1];
rz(-pi) q[2];
rz(-1.6536413) q[3];
sx q[3];
rz(-1.9085064) q[3];
sx q[3];
rz(3.0520671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2485409) q[2];
sx q[2];
rz(-0.38022843) q[2];
sx q[2];
rz(0.60379544) q[2];
rz(2.2327312) q[3];
sx q[3];
rz(-1.7085608) q[3];
sx q[3];
rz(-2.042167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4029082) q[0];
sx q[0];
rz(-1.383902) q[0];
sx q[0];
rz(-0.72147328) q[0];
rz(1.3353222) q[1];
sx q[1];
rz(-1.1386917) q[1];
sx q[1];
rz(-1.8849751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71446705) q[0];
sx q[0];
rz(-0.79475105) q[0];
sx q[0];
rz(0.42596547) q[0];
x q[1];
rz(0.67992731) q[2];
sx q[2];
rz(-0.69831738) q[2];
sx q[2];
rz(0.54496965) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80714204) q[1];
sx q[1];
rz(-0.48105194) q[1];
sx q[1];
rz(-1.0472286) q[1];
x q[2];
rz(-3.0959835) q[3];
sx q[3];
rz(-1.6001587) q[3];
sx q[3];
rz(-1.2330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98715544) q[2];
sx q[2];
rz(-0.64266959) q[2];
sx q[2];
rz(-1.8670234) q[2];
rz(-2.4529723) q[3];
sx q[3];
rz(-1.6504811) q[3];
sx q[3];
rz(-0.0037746599) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0140728) q[0];
sx q[0];
rz(-1.2496244) q[0];
sx q[0];
rz(-2.4697812) q[0];
rz(-1.6067243) q[1];
sx q[1];
rz(-0.22767362) q[1];
sx q[1];
rz(2.6188376) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5005258) q[0];
sx q[0];
rz(-2.4198774) q[0];
sx q[0];
rz(-0.39875491) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7120213) q[2];
sx q[2];
rz(-2.7948973) q[2];
sx q[2];
rz(-0.99682158) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8485496) q[1];
sx q[1];
rz(-0.58401744) q[1];
sx q[1];
rz(2.4021781) q[1];
rz(-pi) q[2];
rz(-1.3626484) q[3];
sx q[3];
rz(-1.2537595) q[3];
sx q[3];
rz(-0.014118346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.25040024) q[2];
sx q[2];
rz(-2.4947391) q[2];
sx q[2];
rz(2.7817173) q[2];
rz(-1.322809) q[3];
sx q[3];
rz(-1.6836932) q[3];
sx q[3];
rz(-3.0944284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5322402) q[0];
sx q[0];
rz(-0.9541963) q[0];
sx q[0];
rz(-0.94616079) q[0];
rz(0.040945176) q[1];
sx q[1];
rz(-1.2766726) q[1];
sx q[1];
rz(-1.2650222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32308043) q[0];
sx q[0];
rz(-1.4671333) q[0];
sx q[0];
rz(2.6635499) q[0];
rz(-pi) q[1];
rz(-1.8184471) q[2];
sx q[2];
rz(-2.0132092) q[2];
sx q[2];
rz(-0.72911791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.62502669) q[1];
sx q[1];
rz(-2.0989162) q[1];
sx q[1];
rz(-2.3305446) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.944421) q[3];
sx q[3];
rz(-1.5253673) q[3];
sx q[3];
rz(1.3988972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1653183) q[2];
sx q[2];
rz(-1.1555305) q[2];
sx q[2];
rz(-0.18132845) q[2];
rz(-0.85339648) q[3];
sx q[3];
rz(-0.80297628) q[3];
sx q[3];
rz(1.8027423) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637909) q[0];
sx q[0];
rz(-1.811857) q[0];
sx q[0];
rz(-1.3743251) q[0];
rz(1.6865267) q[1];
sx q[1];
rz(-1.68119) q[1];
sx q[1];
rz(2.841058) q[1];
rz(2.4886738) q[2];
sx q[2];
rz(-2.3873285) q[2];
sx q[2];
rz(-1.8079226) q[2];
rz(0.28742049) q[3];
sx q[3];
rz(-1.3232373) q[3];
sx q[3];
rz(-0.19178615) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
