// Tiny hash router: #/<tab>/<rest>
function parse() {
  const [tab = 'overview', ...rest] = location.hash.replace(/^#\/?/, '').split('/');
  return { tab: tab || 'overview', rest: rest.map(decodeURIComponent) };
}

export const route = $state(parse());

window.addEventListener('hashchange', () => Object.assign(route, parse()));

export const go = (tab: string, ...rest: string[]) => {
  location.hash = '#/' + [tab, ...rest.map(encodeURIComponent)].join('/');
};
