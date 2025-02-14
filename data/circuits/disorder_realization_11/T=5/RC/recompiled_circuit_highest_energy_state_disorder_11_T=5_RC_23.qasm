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
rz(1.4229245) q[0];
sx q[0];
rz(-2.0473502) q[0];
sx q[0];
rz(-0.25804582) q[0];
rz(2.1482422) q[1];
sx q[1];
rz(-1.300783) q[1];
sx q[1];
rz(-2.8859477) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50798049) q[0];
sx q[0];
rz(-2.2574212) q[0];
sx q[0];
rz(2.344374) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5077816) q[2];
sx q[2];
rz(-1.1828701) q[2];
sx q[2];
rz(2.9390719) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62754831) q[1];
sx q[1];
rz(-1.5020292) q[1];
sx q[1];
rz(-1.3876983) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7840476) q[3];
sx q[3];
rz(-1.0825601) q[3];
sx q[3];
rz(-2.8247716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51097441) q[2];
sx q[2];
rz(-1.7642085) q[2];
sx q[2];
rz(-1.1154491) q[2];
rz(0.014178064) q[3];
sx q[3];
rz(-1.3445798) q[3];
sx q[3];
rz(-2.7052963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2071335) q[0];
sx q[0];
rz(-0.62196982) q[0];
sx q[0];
rz(1.2472664) q[0];
rz(-2.6994052) q[1];
sx q[1];
rz(-1.3860044) q[1];
sx q[1];
rz(2.3242548) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38430518) q[0];
sx q[0];
rz(-1.7624859) q[0];
sx q[0];
rz(2.1038281) q[0];
rz(-pi) q[1];
rz(2.1139718) q[2];
sx q[2];
rz(-0.50293844) q[2];
sx q[2];
rz(0.78928141) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.081233844) q[1];
sx q[1];
rz(-0.66141719) q[1];
sx q[1];
rz(-0.27137406) q[1];
rz(-pi) q[2];
rz(2.3807008) q[3];
sx q[3];
rz(-2.0756654) q[3];
sx q[3];
rz(-2.1438716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3680129) q[2];
sx q[2];
rz(-2.7648338) q[2];
sx q[2];
rz(0.22029857) q[2];
rz(1.9593272) q[3];
sx q[3];
rz(-1.1743098) q[3];
sx q[3];
rz(-3.0784472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0910864) q[0];
sx q[0];
rz(-1.3171221) q[0];
sx q[0];
rz(1.1385981) q[0];
rz(-2.7298722) q[1];
sx q[1];
rz(-0.76985923) q[1];
sx q[1];
rz(2.128111) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6665812) q[0];
sx q[0];
rz(-1.555976) q[0];
sx q[0];
rz(2.8885452) q[0];
rz(-1.1297497) q[2];
sx q[2];
rz(-1.4510321) q[2];
sx q[2];
rz(-2.3074367) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.895993) q[1];
sx q[1];
rz(-2.2618544) q[1];
sx q[1];
rz(-0.57180239) q[1];
rz(-pi) q[2];
rz(-0.13000906) q[3];
sx q[3];
rz(-2.293236) q[3];
sx q[3];
rz(0.60628451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98075214) q[2];
sx q[2];
rz(-1.9858805) q[2];
sx q[2];
rz(1.2551003) q[2];
rz(2.8694425) q[3];
sx q[3];
rz(-0.77747074) q[3];
sx q[3];
rz(3.0009771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.960152) q[0];
sx q[0];
rz(-2.20521) q[0];
sx q[0];
rz(-0.61035672) q[0];
rz(-2.4329674) q[1];
sx q[1];
rz(-2.1253864) q[1];
sx q[1];
rz(-1.5707387) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8886531) q[0];
sx q[0];
rz(-1.4331237) q[0];
sx q[0];
rz(0.15451365) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28684692) q[2];
sx q[2];
rz(-1.45668) q[2];
sx q[2];
rz(-2.508004) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29250568) q[1];
sx q[1];
rz(-1.5226814) q[1];
sx q[1];
rz(2.8699257) q[1];
x q[2];
rz(1.9746154) q[3];
sx q[3];
rz(-2.32956) q[3];
sx q[3];
rz(2.5126575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2107971) q[2];
sx q[2];
rz(-1.0127298) q[2];
sx q[2];
rz(-2.5861758) q[2];
rz(2.4173315) q[3];
sx q[3];
rz(-1.1490425) q[3];
sx q[3];
rz(-2.485062) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.637218) q[0];
sx q[0];
rz(-1.1906304) q[0];
sx q[0];
rz(2.7727238) q[0];
rz(2.8952307) q[1];
sx q[1];
rz(-1.8135704) q[1];
sx q[1];
rz(1.4341644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1911538) q[0];
sx q[0];
rz(-2.3847347) q[0];
sx q[0];
rz(3.1340772) q[0];
rz(-1.7463914) q[2];
sx q[2];
rz(-1.1313952) q[2];
sx q[2];
rz(1.1281769) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1885438) q[1];
sx q[1];
rz(-1.5750139) q[1];
sx q[1];
rz(-0.40178816) q[1];
rz(-3.0415972) q[3];
sx q[3];
rz(-1.7971562) q[3];
sx q[3];
rz(-2.5842427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4904334) q[2];
sx q[2];
rz(-1.8305402) q[2];
sx q[2];
rz(-2.0595713) q[2];
rz(-0.79484445) q[3];
sx q[3];
rz(-1.4719897) q[3];
sx q[3];
rz(1.2580416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.865888) q[0];
sx q[0];
rz(-0.10144932) q[0];
sx q[0];
rz(-0.24359447) q[0];
rz(2.1547735) q[1];
sx q[1];
rz(-0.74816626) q[1];
sx q[1];
rz(2.2056244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2631131) q[0];
sx q[0];
rz(-1.2498444) q[0];
sx q[0];
rz(1.9892938) q[0];
rz(-0.34910874) q[2];
sx q[2];
rz(-0.85452467) q[2];
sx q[2];
rz(-0.30308613) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0985581) q[1];
sx q[1];
rz(-1.1698616) q[1];
sx q[1];
rz(0.75717302) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6020847) q[3];
sx q[3];
rz(-0.60006053) q[3];
sx q[3];
rz(2.4158583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4794856) q[2];
sx q[2];
rz(-1.138849) q[2];
sx q[2];
rz(-2.7916419) q[2];
rz(-0.55772603) q[3];
sx q[3];
rz(-1.2099268) q[3];
sx q[3];
rz(-2.2902655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6996985) q[0];
sx q[0];
rz(-2.5040369) q[0];
sx q[0];
rz(2.8398474) q[0];
rz(0.31983495) q[1];
sx q[1];
rz(-1.539307) q[1];
sx q[1];
rz(1.1357657) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0253191) q[0];
sx q[0];
rz(-2.2145382) q[0];
sx q[0];
rz(-1.5035065) q[0];
rz(-pi) q[1];
rz(2.2280424) q[2];
sx q[2];
rz(-2.3956809) q[2];
sx q[2];
rz(-0.28829703) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7702076) q[1];
sx q[1];
rz(-1.6504297) q[1];
sx q[1];
rz(-2.1712028) q[1];
x q[2];
rz(2.3507422) q[3];
sx q[3];
rz(-1.491427) q[3];
sx q[3];
rz(0.96109238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4096421) q[2];
sx q[2];
rz(-0.66764098) q[2];
sx q[2];
rz(0.14275924) q[2];
rz(-1.7536633) q[3];
sx q[3];
rz(-1.9526491) q[3];
sx q[3];
rz(0.68283844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9614354) q[0];
sx q[0];
rz(-3.1249983) q[0];
sx q[0];
rz(0.55602443) q[0];
rz(0.22932886) q[1];
sx q[1];
rz(-1.2622958) q[1];
sx q[1];
rz(1.3831327) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62263238) q[0];
sx q[0];
rz(-2.3078902) q[0];
sx q[0];
rz(1.6181158) q[0];
rz(-0.45852335) q[2];
sx q[2];
rz(-1.6237729) q[2];
sx q[2];
rz(-0.76901877) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8236716) q[1];
sx q[1];
rz(-1.2814199) q[1];
sx q[1];
rz(-1.9993062) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6881668) q[3];
sx q[3];
rz(-1.9186221) q[3];
sx q[3];
rz(-2.3330101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6250299) q[2];
sx q[2];
rz(-2.8690858) q[2];
sx q[2];
rz(-1.2410835) q[2];
rz(0.062601335) q[3];
sx q[3];
rz(-1.4227941) q[3];
sx q[3];
rz(0.82714287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.1271707) q[0];
sx q[0];
rz(-1.7272471) q[0];
sx q[0];
rz(-0.14778368) q[0];
rz(0.5087018) q[1];
sx q[1];
rz(-1.769442) q[1];
sx q[1];
rz(0.37857372) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0914144) q[0];
sx q[0];
rz(-1.1568406) q[0];
sx q[0];
rz(1.8267836) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0835593) q[2];
sx q[2];
rz(-0.65750018) q[2];
sx q[2];
rz(-3.0160144) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7853773) q[1];
sx q[1];
rz(-1.9058203) q[1];
sx q[1];
rz(-0.32932333) q[1];
x q[2];
rz(1.6949953) q[3];
sx q[3];
rz(-2.327965) q[3];
sx q[3];
rz(0.99909335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.96616894) q[2];
sx q[2];
rz(-2.1899026) q[2];
sx q[2];
rz(1.8625205) q[2];
rz(-1.7133948) q[3];
sx q[3];
rz(-2.1396075) q[3];
sx q[3];
rz(-2.9960347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7487504) q[0];
sx q[0];
rz(-2.5635283) q[0];
sx q[0];
rz(-3.0774935) q[0];
rz(-0.52866689) q[1];
sx q[1];
rz(-1.268498) q[1];
sx q[1];
rz(2.166523) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9720358) q[0];
sx q[0];
rz(-2.0590879) q[0];
sx q[0];
rz(1.8711062) q[0];
rz(-0.56985241) q[2];
sx q[2];
rz(-1.0625216) q[2];
sx q[2];
rz(3.0076671) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84555039) q[1];
sx q[1];
rz(-2.6471889) q[1];
sx q[1];
rz(-2.666996) q[1];
x q[2];
rz(-2.3894839) q[3];
sx q[3];
rz(-0.32882133) q[3];
sx q[3];
rz(3.0260835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.14572445) q[2];
sx q[2];
rz(-0.66103649) q[2];
sx q[2];
rz(0.60689849) q[2];
rz(-0.25035826) q[3];
sx q[3];
rz(-1.4162049) q[3];
sx q[3];
rz(2.6355766) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.09457) q[0];
sx q[0];
rz(-1.9245514) q[0];
sx q[0];
rz(-1.2522329) q[0];
rz(2.6429214) q[1];
sx q[1];
rz(-1.5892727) q[1];
sx q[1];
rz(-1.7472063) q[1];
rz(-0.5947391) q[2];
sx q[2];
rz(-2.1872216) q[2];
sx q[2];
rz(-2.8333153) q[2];
rz(-1.7767033) q[3];
sx q[3];
rz(-1.4510703) q[3];
sx q[3];
rz(2.5616796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
