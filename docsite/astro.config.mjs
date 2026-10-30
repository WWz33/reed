// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// The site is served from https://wwz33.github.io/reed/ — the GitHub Pages
// project path is the Astro `base`. Starlight prefixes its own links with it;
// content links are written relative so they resolve inside either language.
export default defineConfig({
	site: 'https://wwz33.github.io',
	base: '/reed',
	integrations: [
		starlight({
			title: 'reed',
			description: 'Population-genetics statistics and hard filtering, in pure Rust',
			defaultLocale: 'root',
			locales: {
				root: { label: 'English', lang: 'en' },
				'zh-cn': { label: '简体中文', lang: 'zh-CN' },
			},
			social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/WWz33/reed' }],
			sidebar: [
				{
					label: 'Get started',
					translations: { 'zh-CN': '开始' },
					items: ['install', 'quickstart', 'commands'],
				},
				{
					label: 'Command reference',
					translations: { 'zh-CN': '命令参考' },
					items: ['filter', 'stats', 'popg', 'pca', 'distance', 'convert', 'rename', 'index-cmd'],
				},
				{
					label: 'Reference',
					translations: { 'zh-CN': '参考' },
					items: ['expressions', 'outputs', 'library', 'citation'],
				},
			],
		}),
	],
});
